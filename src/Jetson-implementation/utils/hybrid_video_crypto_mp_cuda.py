"""
Hybrid Multi-Processing + CUDA Encryption
Uses existing MP for keystream generation + CUDA for feedback encryption
Best of both worlds: Fast keystream + GPU acceleration
"""

import numpy as np
import sys
import os

# Import the working multi-processing encryption
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


# CUDA kernel for feedback encryption (sequential processing)
CUDA_FEEDBACK_KERNEL = """
__global__ void feedback_encrypt_kernel(
    const unsigned char* plaintext,
    const unsigned char* keystream,
    unsigned char* ciphertext,
    int n
)
{
    // Feedback encryption is inherently sequential
    // Only thread 0 should execute this
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        unsigned char prev = 0;
        for (int i = 0; i < n; i++) {
            ciphertext[i] = (plaintext[i] ^ keystream[i] ^ prev) & 0xFF;
            prev = ciphertext[i];
        }
    }
}

__global__ void feedback_decrypt_kernel(
    const unsigned char* ciphertext,
    const unsigned char* keystream,
    unsigned char* plaintext,
    int n
)
{
    // Feedback decryption is inherently sequential
    // Only thread 0 should execute this
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        unsigned char prev = 0;
        for (int i = 0; i < n; i++) {
            plaintext[i] = (ciphertext[i] ^ keystream[i] ^ prev) & 0xFF;
            prev = ciphertext[i];
        }
    }
}
"""


class HybridVideoEncryptionMPCUDA(HybridVideoEncryptionMP):
    """
    Alias for HybridVideoEncryptionMP - uses parent's optimized CPU implementation.
    
    CUDA acceleration doesn't improve performance for feedback encryption because:
    - Feedback encryption is inherently sequential (each byte depends on previous)
    - CPU with Numba JIT + multiprocessing is faster than sequential GPU execution
    - GPU memory transfers add overhead for sequential operations
    
    This class exists for API compatibility but uses parent's implementation.
    """
    
    def __init__(self, secret_key, frame_width=320, frame_height=240, use_cuda=False):
        # Initialize parent's optimized MP implementation
        super().__init__(frame_width=frame_width, frame_height=frame_height, secret_key=secret_key)
        
        if PYCUDA_AVAILABLE:
            print(f"[INFO] CUDA available: {cuda.Device(0).name()}")
            print(f"[INFO] Using CPU multiprocessing (faster for sequential feedback encryption)")
        else:
            print("[INFO] Using optimized CPU multi-processing")
    
    # Use parent's optimized encrypt_frame() and decrypt_frame() methods
    # No need to override - they use multiprocessing which is faster than sequential CUDA


if __name__ == "__main__":
    import time
    
    print("=" * 70)
    print("Hybrid MP+CUDA Encryption Test")
    print("=" * 70)
    
    # Test 320x240 frame
    test_frame = np.random.randint(0, 256, (240, 320, 3), dtype=np.uint8)
    
    print("\nInitializing MP+CUDA encryption...")
    crypto = HybridVideoEncryptionMPCUDA("test_key", frame_width=320, frame_height=240, use_cuda=True)
    
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
    crypto2 = HybridVideoEncryptionMPCUDA("test_key", frame_width=320, frame_height=240, use_cuda=True)
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
        
        # Test with same instance
        print("\nTesting with same instance...")
        decrypted_same = crypto.decrypt_frame(encrypted)
        if np.array_equal(test_frame, decrypted_same):
            print("✅ Same instance works!")
        else:
            print("❌ Same instance also fails!")
            diff2 = np.abs(test_frame.astype(np.int16) - decrypted_same.astype(np.int16))
            print(f"Max difference: {np.max(diff2)}")
            print(f"Mean difference: {np.mean(diff2):.2f}")
