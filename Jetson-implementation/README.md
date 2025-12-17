# 🧪 Video Encryption System - Complete Guide

Advanced video encryption system with GPU acceleration using CTR mode and CUDA on Jetson Nano.

---

## 📝 Overview

The offline video encryption system provides:
- **GPU-Accelerated Encryption**: CTR+CUDA mode for 16× speedup on Jetson Nano
- **Multiple Encryption Modes**: Hybrid CPU, Multi-processing CPU, or GPU-accelerated
- **Key-based Security**: Strong encryption with chaotic maps (Lorenz, Rössler, Hénon, Tent)
- **CTR Mode**: Counter mode for parallel GPU processing
- **Perfect Reconstruction**: Lossless decryption with correct key

### Performance Comparison

| Mode | Jetson Nano @ 320×240 | Best For |
|------|----------------------|----------|
| **CTR+CUDA** | **82-127 FPS** (7-12ms) | ✅ **Production** (fastest) |
| CPU Multi-processing | 7.6 FPS (131ms) | Limited resources |
| CPU Single-threaded | 2.9 FPS (347ms) | Testing only |

**GPU Speedup**: Up to **16.7× faster** than CPU multi-processing!

---

## 🔐 Encrypt Video File

### ✅ Recommended Usage (GPU-Accelerated - Fastest!)
```bash
python3 encrypt_video_file.py --input test_video.mp4 --output encrypted.enc --mode hybrid_ctr_cuda --key "my_secret_key_2025"
```
**Performance**: 82-127 FPS @ 320×240 on Jetson Nano (16× faster than CPU!)

### Alternative: CPU Multi-processing
```bash
python3 encrypt_video_file.py --input test_video.mp4 --output encrypted.enc --key "my_secret_key_2025" --threads
```
**Performance**: ~7.6 FPS @ 320×240

### Available Options
```bash
python3 encrypt_video_file.py \
  --input test_video.mp4           # Input video file (required)
  --output encrypted.enc            # Output encrypted file (required)
  --key "my_secret_key_2025"        # Secret encryption key (REQUIRED)
  --mode hybrid_ctr_cuda            # GPU mode (fastest) or hybrid (CPU)
  --cuda                            # Enable CUDA acceleration (alternative to mode flag)
  --threads                         # CPU multi-processing (7.6 FPS)
  --width 320                       # Frame width (default: 320)
  --height 240                      # Frame height (default: 240)
```

### Encryption Modes Explained

| Mode | Algorithm | Performance | Best For |
|------|-----------|-------------|----------|
| `hybrid_ctr_cuda` | CTR mode + GPU parallel | 82-127 FPS | ✅ **Production (Jetson)** |
| `hybrid --threads` | Feedback + CPU multiprocessing | 7.6 FPS | Limited GPU access |
| `hybrid` | Feedback + Single CPU | 2.9 FPS | Testing only |

### Expected Output (GPU Mode)
```
============================================================
OFFLINE VIDEO ENCRYPTION
============================================================
Input:  test_video.mp4
Output: encrypted.enc
Mode:   hybrid_ctr_cuda
Size:   320×240
Frames: 300
FPS:    30.00

Initializing encryptor...
[INFO] CTR+CUDA mode - Device: NVIDIA Tegra X1
[INFO] Parallel encryption enabled (CTR mode)
Using CTR+CUDA GPU acceleration (parallel mode)
✓ Initialization: 1.40s

Encrypting 300 frames...
  [ 10.0%] Frame 30/300 | Enc: 12.57ms (79.6 FPS) | Overall: 73.4 FPS
  [ 50.0%] Frame 150/300 | Enc: 11.29ms (88.6 FPS) | Overall: 81.9 FPS
  [100.0%] Frame 300/300 | Enc: 11.22ms (89.1 FPS) | Overall: 82.5 FPS

============================================================
ENCRYPTION COMPLETE
============================================================
Total frames:      300
Total time:        4.04s
Average FPS:       74.34
Avg encrypt time:  11.22ms
File size:         65.93 MB
Save time:         0.40s
Output file:       encrypted.enc
============================================================

✓ Video encrypted successfully!
```

### ⚠️ Security Warning
If you don't provide a `--key`, you'll see:
```
⚠️  WARNING: No encryption key provided!
   Using default key - encryption is NOT secure.
   Use --key 'your_secret_key' for secure encryption.
```

**Always use a custom key for production!**

---

## 🔓 Decrypt Video File

### ✅ Recommended Usage (GPU-Accelerated - Fastest!)
```bash
python3 decrypt_video_file.py --input encrypted.enc --output decrypted.mp4 --mode hybrid_ctr_cuda --key "my_secret_key_2025"
```
**Performance**: 82-127 FPS @ 320×240 on Jetson Nano

### Alternative: CPU Multi-processing
```bash
python3 decrypt_video_file.py --input encrypted.enc --output decrypted.mp4 --key "my_secret_key_2025" --threads
```
**Performance**: ~7.6 FPS @ 320×240

### Available Options
```bash
python3 decrypt_video_file.py \
  --input encrypted.enc             # Input encrypted file (required)
  --output decrypted.mp4            # Output video file (optional - for saving)
  --key "my_secret_key_2025"        # Secret decryption key (MUST match encryption key)
  --mode hybrid_ctr_cuda            # GPU mode (fastest) or hybrid (CPU)
  --cuda                            # Enable CUDA acceleration (alternative to mode flag)
  --threads                         # CPU multi-processing (7.6 FPS)
  --no-display                      # Don't show video while decrypting
```

### Expected Output
```
============================================================
OFFLINE VIDEO DECRYPTION
============================================================
Input:   encrypted.enc
Output:  decrypted.mp4
Display: Yes
Key:     my_secre... (secured)

Loading encrypted data...
✓ Loaded: 4.18 MB in 0.08s
Mode:     hybrid
Size:     320×240
FPS:      30.00
Frames:   150

Initializing decryptor...
Using multi-processing decryption (3 processes)
✓ Initialization: 2.13s

Decrypting 150 frames...
Progress: 100% |████████████████████| 150/150 [00:19<00:00, 7.61 FPS]

============================================================
DECRYPTION COMPLETE
============================================================
Total frames:      150
Total time:        19.71s
Average FPS:       7.61
Avg decrypt time:  131.40ms
Output size:       3.95 MB
Output file:       decrypted.mp4
============================================================
```

### ⚠️ Key Mismatch
Using the wrong key will produce garbage output:
```bash
# Encrypt with one key
python encrypt_video_file.py --input test.mp4 --output test.enc --key "correct_key"

# Decrypt with wrong key = noise/garbage
python decrypt_video_file.py --input test.enc --output wrong.mp4 --key "wrong_key"
# Output: Corrupted video (noise)

# Decrypt with correct key = perfect restoration
python decrypt_video_file.py --input test.enc --output restored.mp4 --key "correct_key"
# Output: Identical to original
```

---

## 📁 Output Files

### 1. Encrypted File (.enc)
Binary file containing:
- Encrypted frame data (pickled Python object)
- Video metadata (resolution, FPS, frame count)
- Cannot be played directly (requires decryption)

**Structure:**
```python
{
    'metadata': {
        'mode': 'hybrid',
        'width': 320,
        'height': 240,
        'fps': 30.0,
        'total_frames': 150,
        'original_file': 'test_video.mp4'
    },
    'frames': [encrypted_frame_1, encrypted_frame_2, ...]
}
```

### 2. Decrypted Video (.mp4)
Standard MP4 video file:
- **Size**: Similar to original
- **Quality**: Lossless (perfect reconstruction with correct key)
- **Playable**: Yes (VLC, mpv, any video player)

---

## 🔬 Technical Details

### CTR Mode vs Feedback Mode

**Why CTR Mode is 16× Faster:**

| Aspect | Feedback Mode (Old) | CTR Mode (New) |
|--------|-------------------|----------------|
| **Encryption** | `C[i] = P[i] ⊕ K[i] ⊕ C[i-1]` | `C[i] = P[i] ⊕ K[i]` |
| **Dependencies** | Sequential (each byte needs previous) | Independent (fully parallel) |
| **GPU Threads** | 1 thread (forced sequential) | 256+ threads (parallel) |
| **Jetson Performance** | 12.6 FPS | 82-127 FPS |
| **Speedup** | 1× baseline | 16.7× faster |


<!-- Add after an appropriate section, e.g., after Section 2 or as a new appendix -->

## Appendix: CTR+CUDA vs CPU Multi-Threading Implementation Comparison

| Feature                | CTR+CUDA (GPU)                                      | CPU Multi-Threading (MP)                |
|------------------------|-----------------------------------------------------|-----------------------------------------|
| **Parallelism**        | Massive: 128 CUDA cores (all pixels in parallel)    | 3 CPU processes (R, G, B in parallel)   |
| **Encryption Mode**    | CTR (Counter) mode: `C[i]=P[i]⊕K[i]` (no feedback) | Feedback mode: `C[i]=P[i]⊕K[i]⊕C[i-1]`  |
| **Permutation**        | GPU kernel (parallel shuffle)                       | NumPy sort/indexing (per process)       |
| **Diffusion**          | XOR with chaotic keystream (CTR, parallel)          | XOR with keystream + feedback (serial)  |
| **Rounds**             | 2 (default, configurable)                           | 2 (default, configurable)               |
| **Keystream**          | Chaotic maps + SHA-256, rotated per round           | Chaotic maps + SHA-256, rotated         |
| **Speed (320×240)**    | ~300–325 FPS                                        | ~7–8 FPS                                |
| **Hardware**           | NVIDIA GPU (Jetson Nano) required                   | Any multi-core CPU                      |
| **Dependencies**       | PyCUDA, CUDA toolkit                                | Python multiprocessing, NumPy           |
| **Initialization**     | ~1.4s (kernel compile, context setup)               | ~1s (process pool, keystreams)          |
| **Best Use**           | Real-time, high-throughput, embedded GPU            | Fallback, batch, no-GPU environments    |
| **Scalability**        | Scales with frame size (GPU parallelism)            | Limited by CPU core count               |
| **Algorithm Security** | Strong (CTR + permutation + chaos)                  | Strong (feedback + permutation + chaos) |
| **Decryption**         | Identical process (CTR is symmetric)                | Identical, but feedback must be reversed|
| **Output**             | Encrypted image/video frames                        | Encrypted image/video frames            |


### Security Analysis

Both modes provide equivalent security:
- **Confusion**: Permutation using chaotic maps (Lorenz, Rössler, Hénon, Tent)
- **Diffusion**: XOR with whitened keystream (SHA-256 based)
- **Multiple Rounds**: Default 3 rounds for thorough mixing
- **Key Derivation**: Chaotic systems seeded with secret key

**CTR mode advantage**: Same security, massively faster on GPU!

### Benchmark on Jetson Nano

```bash
python3 benchmark_performance.py
```

**Expected Results @ 320×240:**
```
Resolution    CPU FPS      CTR+CUDA FPS   Speedup      Real-time?
----------------------------------------------------------------------
160×120       25.32        302.45         11.95×       ✅ Yes
320×240       7.63         127.61         16.72×       ✅ Yes
640×480       2.08         51.23          24.63×       ✅ Yes
```

---

## 📈 Performance Metrics (Jetson Nano)

### Benchmark Results @ 320×240

| Mode | FPS | Time/Frame | GPU Speedup | Use Case |
|------|-----|------------|-------------|----------|
| **CTR+CUDA (GPU)** | **82-127** | **7.8-12ms** | **16.7×** | ✅ **Production** |
| CPU Multi-processing | 7.6 | 131ms | 1× (baseline) | Limited GPU |
| CPU Single-threaded | 2.9 | 347ms | 0.38× | Testing only |

**Key Achievements**:
- 🚀 **127 FPS peak performance** with CTR mode on GPU
- ⚡ **16.7× speedup** over CPU multiprocessing
- 🎯 **Real-time encryption** at 30+ FPS achieved

### Why CTR+CUDA is Faster

| Feature | Feedback Mode (Old) | CTR Mode (New) |
|---------|-------------------|----------------|
| **Parallelization** | ❌ Sequential (byte-by-byte chain) | ✅ Fully parallel |
| **GPU Utilization** | 1 thread only | 256+ threads |
| **Performance** | 12.6 FPS (sequential CUDA) | 127 FPS (parallel CUDA) |
| **Algorithm** | `C[i] = P[i] ⊕ K[i] ⊕ C[i-1]` | `C[i] = P[i] ⊕ K[i]` |

### Resolution Impact

| Resolution | CTR+CUDA | CPU Multi-processing | Speedup |
|------------|----------|---------------------|---------|
| 160×120 | ~300 FPS | ~25 FPS | 12× |
| **320×240** | **82-127 FPS** | **7.6 FPS** | **16.7×** |
| 640×480 | ~50 FPS | ~2.1 FPS | 24× |

**Recommendation**: Use **320×240 with CTR+CUDA** for optimal performance on Jetson Nano.

---

## 🚀 Setup CUDA Environment (First Time Only)

If you get `nvcc: command not found` error:

```bash
# Run this once on Jetson
source fix_cuda_path.sh

# Or make it permanent
./setup_cuda_env.sh
source ~/.bashrc
```

This adds CUDA compiler to your PATH so PyCUDA can compile kernels.

---

## 🐛 Common Issues & Solutions

### Issue 1: "nvcc: No such file or directory" (CUDA)
**Problem**: CUDA compiler not in PATH

**Solution**:
```bash
# Quick fix (current session)
source fix_cuda_path.sh

# Permanent fix
./setup_cuda_env.sh
source ~/.bashrc

# Then run encryption/decryption again
python3 encrypt_video_file.py --input test.mp4 --output test.enc --mode hybrid_ctr_cuda --key "my_key"
```

---

### Issue 2: "ModuleNotFoundError: No module named 'utils'"
**Problem**: Running from wrong directory

**Solution**:
```bash
cd python_implementation/video_encryption
python3 encrypt_video_file.py --input test.mp4 --output test.enc --key "my_key"
```

---

### Issue 3: Low FPS (< 10 FPS)
**Problem**: Not using GPU acceleration

**Solution**:
```bash
# Use CTR+CUDA mode for 16× speedup
python3 encrypt_video_file.py --input test.mp4 --output test.enc --mode hybrid_ctr_cuda --key "my_key"
                                                                   ^^^^^^^^^^^^^^^^^^^^
```

---

### Issue 4: Wrong key / Garbage output
**Problem**: Decryption key doesn't match encryption key

**Example**:
```bash
# Encrypt with key "abc123"
python3 encrypt_video_file.py --input test.mp4 --key "abc123" --mode hybrid_ctr_cuda --output test.enc

# Decrypt with SAME key "abc123" ✅
python3 decrypt_video_file.py --input test.enc --key "abc123" --mode hybrid_ctr_cuda --output restored.mp4

# Decrypt with DIFFERENT key "xyz789" ❌ (produces noise)
python3 decrypt_video_file.py --input test.enc --key "xyz789" --mode hybrid_ctr_cuda --output wrong.mp4
```

**Solution**: Always use the exact same key and mode for both encryption and decryption.

---

### Issue 4: No encryption key warning
**Problem**: Forgot to provide `--key` argument

**Warning**:
```
⚠️  WARNING: No encryption key provided!
   Using default key - encryption is NOT secure.
```

**Solution**:
```bash
# Always provide a custom key
python encrypt_video_file.py --input test.mp4 --output test.enc --key "my_secure_key_2025"
```

---

### Issue 5: File not found error
**Problem**: Input file path is incorrect

**Solution**:
```bash
# Use absolute path or correct relative path
python encrypt_video_file.py --input /full/path/to/video.mp4 --output test.enc --key "my_key"

# Or navigate to correct directory first
cd /path/to/videos
python /path/to/encrypt_video_file.py --input video.mp4 --output test.enc --key "my_key"
```

---

## 🎯 Quick Reference

### Command Templates (GPU Mode - Recommended)

**Encryption (fastest - 82-127 FPS)**:
```bash
python3 encrypt_video_file.py --input VIDEO.mp4 --output ENCRYPTED.enc --mode hybrid_ctr_cuda --key "YOUR_KEY"
```

**Decryption (fastest - 82-127 FPS)**:
```bash
python3 decrypt_video_file.py --input ENCRYPTED.enc --output DECRYPTED.mp4 --mode hybrid_ctr_cuda --key "YOUR_KEY"
```

### Complete Workflow Example

```bash
# 1. Setup CUDA environment (first time only)
source fix_cuda_path.sh

# 2. Encrypt video with GPU
python3 encrypt_video_file.py \
  --input test_video.mp4 \
  --output encrypted.enc \
  --mode hybrid_ctr_cuda \
  --key "my_secure_password_2025"

# Output: 82-127 FPS @ 320×240

# 3. Decrypt video with GPU
python3 decrypt_video_file.py \
  --input encrypted.enc \
  --output decrypted.mp4 \
  --mode hybrid_ctr_cuda \
  --key "my_secure_password_2025"

# Output: Perfect reconstruction at 82-127 FPS
```

### Best Practices
✅ **Use CTR+CUDA mode** for 16× speedup on Jetson Nano  
✅ **Always use a custom key** (not default)  
✅ **Match mode and key** for encryption and decryption  
✅ **Test with short videos first** (5-10 seconds)  
✅ **Use 320×240 resolution** for optimal GPU performance  
✅ **Keep your encryption key secret and secure**  
✅ **Run setup_cuda_env.sh once** to configure CUDA paths  

### Don't
❌ Forget the `--key` argument (insecure without it)  
❌ Use different keys or modes for encryption and decryption  
❌ Share encrypted files without securely sharing the key  
❌ Skip CUDA environment setup (causes nvcc errors)  
❌ Use CPU mode when GPU is available (16× slower)  

---

## 📚 Additional Resources

For more details, see:
- **BENCHMARKS.md** - Complete performance analysis and proof
- **QUICKSTART.md** - Quick start guide with examples
- **OPTIMIZATION_PROOF.md** - Multi-processing optimization details
