# 🧪 Testing Guide

Quick guide for encrypting and decrypting video files using `encrypt_video_file.py` and `decrypt_video_file.py`.

---

## 📝 Overview

The offline video encryption system provides:
- **File-to-file encryption**: Encrypt entire video files and save to disk
- **Key-based security**: Use secret keys for secure encryption (REQUIRED)
- **Multi-processing**: 2.64× speedup with parallel processing
- **Perfect reconstruction**: Lossless decryption with correct key

---

## 🔐 Encrypt Video File

### Basic Usage (Single-threaded)
```bash
python encrypt_video_file.py --input test_video.mp4 --output encrypted.enc --key "my_secret_key_2025"
```

### Recommended Usage (Multi-processing)
```bash
python encrypt_video_file.py --input test_video.mp4 --output encrypted.enc --key "my_secret_key_2025" --threads
```

### Available Options
```bash
python encrypt_video_file.py \
  --input test_video.mp4           # Input video file (required)
  --output encrypted.enc            # Output encrypted file (required)
  --key "my_secret_key_2025"        # Secret encryption key (STRONGLY RECOMMENDED)
  --threads                         # Enable multi-processing (2.64x faster)
  --width 320                       # Frame width (default: 320)
  --height 240                      # Frame height (default: 240)
  --mode hybrid                     # Encryption mode (hybrid/fast/lightweight)
```

### Expected Output
```
============================================================
OFFLINE VIDEO ENCRYPTION
============================================================
Input:  test_video.mp4
Output: encrypted.enc
Mode:   hybrid_MT
Size:   320×240
Key:    my_secre... (secured)
Threads: Multi-threaded (3 workers for RGB)
Frames: 150
FPS:    30.00

Initializing encryptor...
Using multi-processing (3 processes) for true parallelism
✓ Initialization: 2.15s

Encrypting 150 frames...
Progress: 100% |████████████████████| 150/150 [00:19<00:00, 7.59 FPS]

============================================================
ENCRYPTION COMPLETE
============================================================
Total frames:      150
Total time:        19.76s
Average FPS:       7.59
Avg encrypt time:  131.73ms
File size:         4.18 MB
Save time:         0.12s
Output file:       encrypted.enc
============================================================
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

### Basic Usage (Display only)
```bash
python decrypt_video_file.py --input encrypted.enc --key "my_secret_key_2025"
```

### Recommended Usage (Save to file with multi-processing)
```bash
python decrypt_video_file.py --input encrypted.enc --output decrypted.mp4 --key "my_secret_key_2025" --threads
```

### Available Options
```bash
python decrypt_video_file.py \
  --input encrypted.enc             # Input encrypted file (required)
  --output decrypted.mp4            # Output video file (optional - for saving)
  --key "my_secret_key_2025"        # Secret decryption key (MUST match encryption key)
  --threads                         # Enable multi-processing (2.64x faster)
  --no-display                      # Don't show video while decrypting
  --mode hybrid                     # Decryption mode (must match encryption)
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

## 🔍 Complete Workflow Example

### Step 1: Encrypt with custom key
```bash
python encrypt_video_file.py \
  --input test_video.mp4 \
  --output test.enc \
  --key "my_secure_password_2025" \
  --threads
```

### Step 2: Decrypt with same key
```bash
python decrypt_video_file.py \
  --input test.enc \
  --output restored.mp4 \
  --key "my_secure_password_2025" \
  --threads
```

### Step 3: Verify restoration
```bash
# Play original
vlc test_video.mp4

# Play restored (should be identical)
vlc restored.mp4

# Or compare frame-by-frame (advanced)
ffmpeg -i test_video.mp4 -i restored.mp4 -filter_complex psnr -f null -
```

---

## 📈 Performance Metrics

### Expected Performance (Jetson Nano @ 320×240)

| Method | FPS | Time/Frame | Use Case |
|--------|-----|------------|----------|
| **Single-threaded** | 2.88 | 347ms | Testing only |
| **Multi-processing** | 7.59 | 132ms | ✅ **Recommended** |

**Speedup**: 2.64× faster with `--threads` flag

### Resolution Impact

| Resolution | Single-threaded | Multi-processing | Speedup |
|------------|-----------------|------------------|---------|
| 160×120 | ~10 FPS | ~25 FPS | 2.5× |
| 320×240 | 2.88 FPS | 7.59 FPS | 2.64× |
| 640×480 | ~0.8 FPS | ~2.1 FPS | 2.6× |

**Recommendation**: Use 320×240 for best balance of speed and quality on Jetson Nano.

---

## 🐛 Common Issues & Solutions

### Issue 1: "ModuleNotFoundError: No module named 'utils'"
**Problem**: Running from wrong directory

**Solution**:
```bash
cd python_implementation/video_encryption
python encrypt_video_file.py --input test.mp4 --output test.enc --key "my_key"
```

---

### Issue 2: Low FPS (< 3 FPS)
**Problem**: Not using multi-processing

**Solution**:
```bash
# Add --threads flag
python encrypt_video_file.py --input test.mp4 --output test.enc --key "my_key" --threads
                                                                                 ^^^^^^^^
```

---

### Issue 3: Wrong key / Garbage output
**Problem**: Decryption key doesn't match encryption key

**Example**:
```bash
# Encrypt with key "abc123"
python encrypt_video_file.py --input test.mp4 --key "abc123" --output test.enc

# Decrypt with SAME key "abc123" ✅
python decrypt_video_file.py --input test.enc --key "abc123" --output restored.mp4

# Decrypt with DIFFERENT key "xyz789" ❌ (produces noise)
python decrypt_video_file.py --input test.enc --key "xyz789" --output wrong.mp4
```

**Solution**: Always use the exact same key for both encryption and decryption.

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

### Command Templates

**Encryption (recommended)**:
```bash
python encrypt_video_file.py --input VIDEO.mp4 --output ENCRYPTED.enc --key "YOUR_KEY" --threads
```

**Decryption (recommended)**:
```bash
python decrypt_video_file.py --input ENCRYPTED.enc --output DECRYPTED.mp4 --key "YOUR_KEY" --threads
```

### Best Practices
✅ **Always use a custom key** (not default)  
✅ **Always use `--threads`** for production (2.64× faster)  
✅ **Test with short videos first** (5-10 seconds)  
✅ **Use 320×240 resolution** for Jetson Nano  
✅ **Keep your encryption key secret and secure**  

### Don't
❌ Forget the `--key` argument (insecure without it)  
❌ Use different keys for encryption and decryption  
❌ Share encrypted files without securely sharing the key  
❌ Encrypt large videos without testing performance first  

---

## 📚 Additional Resources

For more details, see:
- **BENCHMARKS.md** - Complete performance analysis and proof
- **QUICKSTART.md** - Quick start guide with examples
- **OPTIMIZATION_PROOF.md** - Multi-processing optimization details
