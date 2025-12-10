"""
Hybrid Video Encryption with CTR Mode + CUDA
Uses Counter (CTR) mode instead of feedback mode for true parallel GPU execution

CTR Mode Benefits:
- Each byte encrypted independently: ciphertext[i] = plaintext[i] XOR keystream[i]
- No sequential dependencies - perfect for GPU parallelization
- Same security with confusion (permutation) + diffusion (XOR with keystream)

Algorithm:
1. Generate chaotic keystreams (reuse from parent)
2. For each round:
   - Apply permutation (confusion)
   - Apply CTR encryption: parallel XOR with rotated keystream (diffusion)
"""

import numpy as np
import sys
import os

# Import the working multi-processing encryption for keystream generation
sys.path.insert(0, os.path.dirname(__file__))
from hybrid_video_crypto_mp import HybridVideoEncryptionMP

# Try to import CUDA
try:
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    PYCUDA_AVAILABLE = True
except ImportError:
    PYCUDA_AVAILABLE = False
    print("[WARNING] PyCUDA not available - falling back to CPU")


# CUDA kernels for parallel CTR mode encryption
CUDA_CTR_KERNEL = """
__global__ void ctr_encrypt_kernel(
    const unsigned char* plaintext,
    const unsigned char* keystream,
    unsigned char* ciphertext,
    int n
)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < n) {
        // CTR mode: simple XOR, no dependencies
        ciphertext[idx] = (plaintext[idx] ^ keystream[idx]) & 0xFF;
    }
}

__global__ void ctr_decrypt_kernel(
    const unsigned char* ciphertext,
    const unsigned char* keystream,
    unsigned char* plaintext,
    int n
)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < n) {
        // CTR mode: decryption is identical to encryption (XOR is reversible)
        plaintext[idx] = (ciphertext[idx] ^ keystream[idx]) & 0xFF;
    }
}

__global__ void apply_permutation_kernel(
    const unsigned char* input,
    unsigned char* output,
    const int* permutation,
    int n
)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < n) {
        output[idx] = input[permutation[idx]];
    }
}

__global__ void apply_inverse_permutation_kernel(
    const unsigned char* input,
    unsigned char* output,
    const int* inv_permutation,
    int n
)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < n) {
        output[idx] = input[inv_permutation[idx]];
    }
}
"""

# Global kernel cache to avoid recompilation
_CUDA_MODULE_CACHE = None


class HybridVideoEncryptionCTRCUDA(HybridVideoEncryptionMP):
    """
    Hybrid encryption using CTR mode + CUDA for true parallel execution.
    
    Key differences from feedback mode:
    - CTR mode: ciphertext[i] = plaintext[i] XOR keystream[i] (parallel)
    - Feedback mode: ciphertext[i] = plaintext[i] XOR keystream[i] XOR ciphertext[i-1] (sequential)
    
    This allows full GPU parallelization while maintaining security through:
    - Multiple rounds of permutation + CTR encryption
    - Chaotic keystream generation (unpredictable)
    - Keystream rotation per round
    """
    
    def __init__(self, secret_key, frame_width=320, frame_height=240, use_cuda=True):
        # Initialize parent to generate chaotic keystreams
        super().__init__(frame_width=frame_width, frame_height=frame_height, secret_key=secret_key)
        
        self.use_cuda = use_cuda and PYCUDA_AVAILABLE
        
        if self.use_cuda:
            self._initialize_cuda()
            print(f"[INFO] CTR+CUDA mode - Device: {cuda.Device(0).name()}")
            print(f"[INFO] Parallel encryption enabled (CTR mode)")
        else:
            print("[INFO] CTR+CPU mode (CUDA not available)")
    
    def _initialize_cuda(self):
        """Compile CUDA kernels (cached to avoid recompilation)"""
        global _CUDA_MODULE_CACHE
        
        if _CUDA_MODULE_CACHE is None:
            try:
                import os
                # Set cache directory to avoid nvcc dependency
                cache_dir = os.path.expanduser("~/.pycuda_cache")
                os.makedirs(cache_dir, exist_ok=True)
                
                # Compile with caching enabled
                _CUDA_MODULE_CACHE = SourceModule(
                    CUDA_CTR_KERNEL,
                    cache_dir=cache_dir,
                    no_extern_c=False
                )
            except Exception as e:
                print(f"[WARNING] CUDA compilation failed: {e}")
                print("[INFO] Falling back to CPU mode")
                self.use_cuda = False
                return
        
        self.cuda_module = _CUDA_MODULE_CACHE
        self.ctr_encrypt_gpu = self.cuda_module.get_function("ctr_encrypt_kernel")
        self.ctr_decrypt_gpu = self.cuda_module.get_function("ctr_decrypt_kernel")
        self.apply_perm_gpu = self.cuda_module.get_function("apply_permutation_kernel")
        self.apply_inv_perm_gpu = self.cuda_module.get_function("apply_inverse_permutation_kernel")
        
        # GPU buffers (allocated on first use)
        self.gpu_buffers = {}
        self.block_size = 256
    
    def _allocate_gpu_buffers(self, size):
        """Allocate reusable GPU memory for data and permutations"""
        if size not in self.gpu_buffers:
            self.gpu_buffers[size] = {
                'input': cuda.mem_alloc(size),
                'output': cuda.mem_alloc(size),
                'temp': cuda.mem_alloc(size),
                'keystream': cuda.mem_alloc(size),
                'perm': cuda.mem_alloc(size * 4),  # int32 = 4 bytes
                'inv_perm': cuda.mem_alloc(size * 4)
            }
    
    def _ctr_encrypt_cuda(self, plaintext, keystream):
        """GPU-accelerated parallel CTR encryption"""
        n = len(plaintext)
        self._allocate_gpu_buffers(n)
        
        d_input = self.gpu_buffers[n]['input']
        d_keystream = self.gpu_buffers[n]['keystream']
        d_output = self.gpu_buffers[n]['output']
        
        # Transfer to GPU
        cuda.memcpy_htod(d_input, plaintext)
        cuda.memcpy_htod(d_keystream, keystream)
        
        # Launch parallel kernel
        grid_size = (n + self.block_size - 1) // self.block_size
        self.ctr_encrypt_gpu(
            d_input, d_keystream, d_output, np.int32(n),
            block=(self.block_size, 1, 1),
            grid=(grid_size, 1)
        )
        
        # Get result
        ciphertext = np.empty(n, dtype=np.uint8)
        cuda.memcpy_dtoh(ciphertext, d_output)
        return ciphertext
    
    def _ctr_decrypt_cuda(self, ciphertext, keystream):
        """GPU-accelerated parallel CTR decryption (identical to encryption)"""
        n = len(ciphertext)
        self._allocate_gpu_buffers(n)
        
        d_input = self.gpu_buffers[n]['input']
        d_keystream = self.gpu_buffers[n]['keystream']
        d_output = self.gpu_buffers[n]['output']
        
        cuda.memcpy_htod(d_input, ciphertext)
        cuda.memcpy_htod(d_keystream, keystream)
        
        # Launch parallel kernel
        grid_size = (n + self.block_size - 1) // self.block_size
        self.ctr_decrypt_gpu(
            d_input, d_keystream, d_output, np.int32(n),
            block=(self.block_size, 1, 1),
            grid=(grid_size, 1)
        )
        
        plaintext = np.empty(n, dtype=np.uint8)
        cuda.memcpy_dtoh(plaintext, d_output)
        return plaintext
    
    def _apply_permutation_cuda(self, data, permutation):
        """GPU-accelerated parallel permutation"""
        n = len(data)
        self._allocate_gpu_buffers(n)
        
        d_input = self.gpu_buffers[n]['input']
        d_output = self.gpu_buffers[n]['output']
        d_perm = self.gpu_buffers[n]['perm']
        
        # Transfer data and permutation
        cuda.memcpy_htod(d_input, data)
        cuda.memcpy_htod(d_perm, permutation.astype(np.int32))
        
        # Launch parallel kernel
        grid_size = (n + self.block_size - 1) // self.block_size
        self.apply_perm_gpu(
            d_input, d_output, d_perm, np.int32(n),
            block=(self.block_size, 1, 1),
            grid=(grid_size, 1)
        )
        
        result = np.empty(n, dtype=np.uint8)
        cuda.memcpy_dtoh(result, d_output)
        return result
    
    def _apply_inverse_permutation_cuda(self, data, inv_permutation):
        """GPU-accelerated parallel inverse permutation"""
        n = len(data)
        self._allocate_gpu_buffers(n)
        
        d_input = self.gpu_buffers[n]['input']
        d_output = self.gpu_buffers[n]['output']
        d_inv_perm = self.gpu_buffers[n]['inv_perm']
        
        # Transfer data and inverse permutation
        cuda.memcpy_htod(d_input, data)
        cuda.memcpy_htod(d_inv_perm, inv_permutation.astype(np.int32))
        
        # Launch parallel kernel
        grid_size = (n + self.block_size - 1) // self.block_size
        self.apply_inv_perm_gpu(
            d_input, d_output, d_inv_perm, np.int32(n),
            block=(self.block_size, 1, 1),
            grid=(grid_size, 1)
        )
        
        result = np.empty(n, dtype=np.uint8)
        cuda.memcpy_dtoh(result, d_output)
        return result
    
    def encrypt_frame(self, frame):
        """
        Encrypt frame using CTR mode + CUDA
        Algorithm: For each round: permutation -> parallel CTR encryption
        """
        h, w = frame.shape[:2]
        
        # Get the RGB channels
        r_channel = frame[:, :, 0].flatten()
        g_channel = frame[:, :, 1].flatten()
        b_channel = frame[:, :, 2].flatten()
        
        # Process each channel
        results = []
        for channel_data, perm, keystream in [
            (r_channel, self.perm_r, self.ks_r),
            (g_channel, self.perm_g, self.ks_g),
            (b_channel, self.perm_b, self.ks_b)
        ]:
            result = channel_data.copy()
            
            for rnd in range(self.num_rounds):
                # Confusion: permutation (parallel on GPU)
                if self.use_cuda:
                    result = self._apply_permutation_cuda(result, perm[:len(result)])
                else:
                    result = result[perm[:len(result)]]
                
                # Diffusion: CTR encryption with keystream rotation (parallel on GPU)
                offset = rnd * 13
                ks = np.roll(keystream, -offset)[:len(result)]
                
                if self.use_cuda:
                    result = self._ctr_encrypt_cuda(result, ks)
                else:
                    # CPU fallback: simple XOR
                    result = (result ^ ks) & 0xFF
            
            results.append(result)
        
        # Reshape back to image
        encrypted_frame = np.stack([
            results[0].reshape(h, w),
            results[1].reshape(h, w),
            results[2].reshape(h, w)
        ], axis=2)
        
        return encrypted_frame
    
    def decrypt_frame(self, encrypted_frame):
        """
        Decrypt frame using CTR mode + CUDA
        Algorithm: For each round (reversed): parallel CTR decryption -> inverse permutation
        """
        h, w = encrypted_frame.shape[:2]
        
        # Get the RGB channels
        r_channel = encrypted_frame[:, :, 0].flatten()
        g_channel = encrypted_frame[:, :, 1].flatten()
        b_channel = encrypted_frame[:, :, 2].flatten()
        
        # Process each channel (reverse order)
        results = []
        for channel_data, inv_perm, keystream in [
            (r_channel, self.inv_perm_r, self.ks_r),
            (g_channel, self.inv_perm_g, self.ks_g),
            (b_channel, self.inv_perm_b, self.ks_b)
        ]:
            result = channel_data.copy()
            
            # Reverse rounds
            for rnd in range(self.num_rounds - 1, -1, -1):
                # Reverse diffusion: CTR decryption (parallel on GPU)
                offset = rnd * 13
                ks = np.roll(keystream, -offset)[:len(result)]
                
                if self.use_cuda:
                    result = self._ctr_decrypt_cuda(result, ks)
                else:
                    # CPU fallback: simple XOR (same as encryption)
                    result = (result ^ ks) & 0xFF
                
                # Reverse confusion: inverse permutation (parallel on GPU)
                if self.use_cuda:
                    result = self._apply_inverse_permutation_cuda(result, inv_perm[:len(result)])
                else:
                    result = result[inv_perm[:len(result)]]
            
            results.append(result)
        
        # Reshape back to image
        decrypted_frame = np.stack([
            results[0].reshape(h, w),
            results[1].reshape(h, w),
            results[2].reshape(h, w)
        ], axis=2)
        
        return decrypted_frame


if __name__ == "__main__":
    import time
    
    print("="*70)
    print("Hybrid CTR+CUDA Encryption Test")
    print("="*70)
    
    # Test 320x240 frame
    test_frame = np.random.randint(0, 256, (240, 320, 3), dtype=np.uint8)
    
    print("\nInitializing CTR+CUDA encryption...")
    crypto = HybridVideoEncryptionCTRCUDA("test_key", frame_width=320, frame_height=240, use_cuda=True)
    
    print("\nWarm-up (5 frames)...")
    for _ in range(5):
        _ = crypto.encrypt_frame(test_frame)
    
    print("\nBenchmarking 30 frames...")
    start = time.time()
    for _ in range(30):
        encrypted = crypto.encrypt_frame(test_frame)
    elapsed = time.time() - start
    
    fps = 30 / elapsed
    ms_per_frame = (elapsed / 30) * 1000
    
    print(f"\n{'='*70}")
    print(f"Results: 320×240")
    print(f"{'='*70}")
    print(f"Total time: {elapsed:.2f}s")
    print(f"FPS: {fps:.2f}")
    print(f"ms/frame: {ms_per_frame:.1f}ms")
    print(f"{'='*70}")
    
    # Verify correctness
    print("\nVerifying correctness...")
    crypto2 = HybridVideoEncryptionCTRCUDA("test_key", frame_width=320, frame_height=240, use_cuda=True)
    decrypted = crypto2.decrypt_frame(encrypted)
    
    if np.array_equal(test_frame, decrypted):
        print("✅ Encryption/decryption verified!")
    else:
        print("❌ Verification failed!")
        # Debug info
        diff = np.abs(test_frame.astype(np.int16) - decrypted.astype(np.int16))
        print(f"Max difference: {np.max(diff)}")
        print(f"Mean difference: {np.mean(diff):.2f}")
        print(f"Pixels different: {np.sum(diff > 0)} / {test_frame.size}")
