#!/usr/bin/env python3
"""
CPU Multi-Threading vs CUDA CTR Mode Performance Test
Direct comparison between the two implementations tested on Jetson Nano
"""

import numpy as np
import time
import sys
import os

# Add utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'utils'))

print("=" * 80)
print("CPU MULTI-THREADING vs CUDA CTR MODE - PERFORMANCE TEST")
print("=" * 80)

# Configuration
# Use a simple key that won't cause numerical overflow in chaotic maps
SECRET_KEY = "mysecretkey"
N_FRAMES = 30
RESOLUTIONS = [
    (160, 120, "160×120"),
    (320, 240, "320×240"),
    (640, 480, "640×480"),
]

# Check implementations
print("\n[1/3] Checking available implementations...")

# CPU Multi-Threading
CPU_AVAILABLE = False
try:
    from hybrid_video_crypto_mp import HybridVideoEncryptionMP
    CPU_AVAILABLE = True
    print("   ✅ CPU Multi-Threading: hybrid_video_crypto_mp.py")
except ImportError as e:
    print(f"   ❌ CPU Multi-Threading not found: {e}")

# CUDA CTR
CUDA_AVAILABLE = False
try:
    from hybrid_video_crypto_ctr_cuda import HybridVideoEncryptionCTRCUDA
    import pycuda.driver as cuda
    import pycuda.autoinit
    CUDA_AVAILABLE = True
    print(f"   ✅ CUDA CTR Mode: hybrid_video_crypto_ctr_cuda.py")
    print(f"      GPU Device: {cuda.Device(0).name()}")
    print(f"      Compute Capability: {cuda.Device(0).compute_capability()}")
except ImportError as e:
    print(f"   ❌ CUDA CTR not available: {e}")

if not CPU_AVAILABLE and not CUDA_AVAILABLE:
    print("\n❌ No implementations available. Exiting.")
    sys.exit(1)

# Results storage
all_results = []

# Run tests
for width, height, name in RESOLUTIONS:
    print("\n" + "=" * 80)
    print(f"[2/3] Testing Resolution: {name} ({height}×{width})")
    print("=" * 80)
    
    # Generate test frame
    test_frame = np.random.randint(0, 256, (height, width, 3), dtype=np.uint8)
    total_pixels = height * width
    
    result = {
        'resolution': name,
        'width': width,
        'height': height,
        'pixels': total_pixels,
        'cpu_mp': None,
        'cuda_ctr': None,
    }
    
    # ========================================================================
    # Test 1: CPU Multi-Threading
    # ========================================================================
    if CPU_AVAILABLE:
        print(f"\n[Test 1/2] CPU Multi-Threading (3 processes)...")
        print("-" * 80)
        
        # Initialize
        crypto_cpu = HybridVideoEncryptionMP(
            frame_width=width,
            frame_height=height,
            secret_key=SECRET_KEY
        )
        
        # Warm-up (3 frames)
        print("   Warm-up: 3 frames...")
        for _ in range(3):
            _ = crypto_cpu.encrypt_frame(test_frame)
        
        # Benchmark
        print(f"   Benchmarking: {N_FRAMES} frames...")
        start_time = time.time()
        
        for i in range(N_FRAMES):
            encrypted = crypto_cpu.encrypt_frame(test_frame)
            if (i + 1) % 10 == 0:
                print(f"      Progress: {i+1}/{N_FRAMES} frames", end='\r')
        
        elapsed = time.time() - start_time
        fps = N_FRAMES / elapsed
        ms_per_frame = (elapsed / N_FRAMES) * 1000
        
        result['cpu_mp'] = {
            'time': elapsed,
            'fps': fps,
            'ms_per_frame': ms_per_frame
        }
        
        print(f"\n   ✅ Results:")
        print(f"      Total time:  {elapsed:.2f} seconds")
        print(f"      FPS:         {fps:.2f}")
        print(f"      ms/frame:    {ms_per_frame:.2f}ms")
        
        # Save encrypted frame for verification
        cpu_encrypted = encrypted.copy()
    
    # ========================================================================
    # Test 2: CUDA CTR Mode
    # ========================================================================
    if CUDA_AVAILABLE:
        print(f"\n[Test 2/2] CUDA CTR Mode (GPU parallel)...")
        print("-" * 80)
        
        # Initialize
        crypto_cuda = HybridVideoEncryptionCTRCUDA(
            secret_key=SECRET_KEY,
            frame_width=width,
            frame_height=height,
            use_cuda=True
        )
        
        # Warm-up (5 frames for GPU)
        print("   Warm-up: 5 frames (GPU kernel compilation)...")
        for _ in range(5):
            _ = crypto_cuda.encrypt_frame(test_frame)
        
        # Benchmark
        print(f"   Benchmarking: {N_FRAMES} frames...")
        start_time = time.time()
        
        for i in range(N_FRAMES):
            encrypted = crypto_cuda.encrypt_frame(test_frame)
            if (i + 1) % 10 == 0:
                print(f"      Progress: {i+1}/{N_FRAMES} frames", end='\r')
        
        elapsed = time.time() - start_time
        fps = N_FRAMES / elapsed
        ms_per_frame = (elapsed / N_FRAMES) * 1000
        
        result['cuda_ctr'] = {
            'time': elapsed,
            'fps': fps,
            'ms_per_frame': ms_per_frame
        }
        
        print(f"\n   ✅ Results:")
        print(f"      Total time:  {elapsed:.2f} seconds")
        print(f"      FPS:         {fps:.2f}")
        print(f"      ms/frame:    {ms_per_frame:.2f}ms")
        
        # Calculate speedup
        if result['cpu_mp']:
            speedup = result['cuda_ctr']['fps'] / result['cpu_mp']['fps']
            print(f"\n   🚀 GPU Speedup: {speedup:.2f}× faster than CPU multi-threading")
        
        # Save encrypted frame
        cuda_encrypted = encrypted.copy()
    
    # ========================================================================
    # Verification
    # ========================================================================
    if CPU_AVAILABLE and CUDA_AVAILABLE:
        print(f"\n[3/3] Verification...")
        print("-" * 80)
        
        # Both implementations should produce different encrypted outputs
        # (different encryption modes), but both should decrypt correctly
        
        # Decrypt CPU encrypted frame
        cpu_decrypted = crypto_cpu.decrypt_frame(cpu_encrypted)
        cpu_match = np.array_equal(test_frame, cpu_decrypted)
        
        # Decrypt CUDA encrypted frame
        cuda_decrypted = crypto_cuda.decrypt_frame(cuda_encrypted)
        cuda_match = np.array_equal(test_frame, cuda_decrypted)
        
        print(f"   CPU encryption/decryption:  {'✅ PASS' if cpu_match else '❌ FAIL'}")
        print(f"   CUDA encryption/decryption: {'✅ PASS' if cuda_match else '❌ FAIL'}")
        
        if cpu_match and cuda_match:
            print("   ✅ Both implementations verified successfully!")
        else:
            print("   ⚠️  Verification issue detected")
            if not cpu_match:
                diff = np.abs(test_frame.astype(np.int16) - cpu_decrypted.astype(np.int16))
                print(f"      CPU: {np.sum(diff > 0)} pixels differ, max diff: {np.max(diff)}")
            if not cuda_match:
                diff = np.abs(test_frame.astype(np.int16) - cuda_decrypted.astype(np.int16))
                print(f"      CUDA: {np.sum(diff > 0)} pixels differ, max diff: {np.max(diff)}")
    
    all_results.append(result)

# ============================================================================
# Summary Table
# ============================================================================
print("\n" + "=" * 80)
print("PERFORMANCE SUMMARY")
print("=" * 80)

# Header
print(f"\n{'Resolution':<12} {'Pixels':<10} {'CPU (MP)':<14} {'CUDA (CTR)':<14} {'Speedup':<12} {'30 FPS?'}")
print("-" * 80)

# Rows
for r in all_results:
    res = r['resolution']
    pixels = f"{r['pixels']:,}"
    
    if r['cpu_mp']:
        cpu_fps = f"{r['cpu_mp']['fps']:.2f} FPS"
        cpu_ms = f"({r['cpu_mp']['ms_per_frame']:.1f}ms)"
    else:
        cpu_fps = "N/A"
        cpu_ms = ""
    
    if r['cuda_ctr']:
        cuda_fps = f"{r['cuda_ctr']['fps']:.2f} FPS"
        cuda_ms = f"({r['cuda_ctr']['ms_per_frame']:.1f}ms)"
    else:
        cuda_fps = "N/A"
        cuda_ms = ""
    
    if r['cpu_mp'] and r['cuda_ctr']:
        speedup = r['cuda_ctr']['fps'] / r['cpu_mp']['fps']
        speedup_str = f"{speedup:.2f}×"
        
        # Check if CUDA can do real-time (30 FPS)
        if r['cuda_ctr']['fps'] >= 30:
            realtime = "✅ Yes"
        elif r['cuda_ctr']['fps'] >= 15:
            realtime = "🟡 Partial"
        else:
            realtime = "❌ No"
    else:
        speedup_str = "N/A"
        realtime = "N/A"
    
    print(f"{res:<12} {pixels:<10} {cpu_fps:<14} {cuda_fps:<14} {speedup_str:<12} {realtime}")

print("-" * 80)

# ============================================================================
# Detailed Breakdown (320×240)
# ============================================================================
if len(all_results) >= 2:
    print("\n" + "=" * 80)
    print("DETAILED ANALYSIS - 320×240 (Standard Resolution)")
    print("=" * 80)
    
    std_result = all_results[1]  # 320×240 is index 1
    
    if std_result['cpu_mp'] and std_result['cuda_ctr']:
        cpu_data = std_result['cpu_mp']
        cuda_data = std_result['cuda_ctr']
        speedup = cuda_data['fps'] / cpu_data['fps']
        
        print(f"\nCPU Multi-Threading (3 processes):")
        print(f"   FPS:          {cpu_data['fps']:.2f}")
        print(f"   ms/frame:     {cpu_data['ms_per_frame']:.2f}ms")
        print(f"   Total time:   {cpu_data['time']:.2f}s for {N_FRAMES} frames")
        
        print(f"\nCUDA CTR Mode (GPU parallel):")
        print(f"   FPS:          {cuda_data['fps']:.2f}")
        print(f"   ms/frame:     {cuda_data['ms_per_frame']:.2f}ms")
        print(f"   Total time:   {cuda_data['time']:.2f}s for {N_FRAMES} frames")
        
        print(f"\nSpeedup Analysis:")
        print(f"   GPU is {speedup:.2f}× faster than CPU multi-threading")
        print(f"   Time saved: {cpu_data['time'] - cuda_data['time']:.2f}s per {N_FRAMES} frames")
        
        # Real-time capability
        cpu_streams_30fps = int(cpu_data['fps'] / 30)
        cuda_streams_30fps = int(cuda_data['fps'] / 30)
        
        print(f"\nReal-time Video Capability (30 FPS target):")
        print(f"   CPU can handle:  {cpu_streams_30fps} simultaneous stream(s)")
        print(f"   GPU can handle:  {cuda_streams_30fps} simultaneous stream(s)")
        print(f"   GPU advantage:   {cuda_streams_30fps - cpu_streams_30fps} more streams")

# ============================================================================
# Recommendations
# ============================================================================
print("\n" + "=" * 80)
print("RECOMMENDATIONS")
print("=" * 80)

if CUDA_AVAILABLE and all_results[1]['cuda_ctr']:
    cuda_320_fps = all_results[1]['cuda_ctr']['fps']
    
    print(f"\n✅ GPU Performance Assessment (320×240):")
    if cuda_320_fps >= 100:
        print(f"   EXCELLENT: {cuda_320_fps:.0f} FPS - production ready!")
        print(f"   → Can handle {int(cuda_320_fps/30)} concurrent 30 FPS video streams")
        print(f"   → Recommended for real-time applications")
    elif cuda_320_fps >= 30:
        print(f"   GOOD: {cuda_320_fps:.0f} FPS - meets real-time threshold")
        print(f"   → Suitable for single-stream real-time encryption")
    else:
        print(f"   MARGINAL: {cuda_320_fps:.0f} FPS - below real-time threshold")
        print(f"   → Consider optimization or lower resolution")
    
    print(f"\n📊 Implementation Choice:")
    if all_results[1]['cpu_mp']:
        cpu_320_fps = all_results[1]['cpu_mp']['fps']
        speedup = cuda_320_fps / cpu_320_fps
        
        if speedup > 20:
            print(f"   → Use CUDA CTR mode (hybrid_video_crypto_ctr_cuda.py)")
            print(f"   → {speedup:.0f}× speedup justifies GPU overhead")
        elif speedup > 5:
            print(f"   → Use CUDA CTR mode for production")
            print(f"   → {speedup:.0f}× speedup provides good benefit")
        else:
            print(f"   → Consider CPU multi-threading for simpler deployment")
            print(f"   → {speedup:.0f}× speedup may not justify GPU complexity")
    else:
        print(f"   → Use CUDA CTR mode (hybrid_video_crypto_ctr_cuda.py)")
else:
    print(f"\n⚠️  GPU testing not available")
    print(f"   → Install PyCUDA to enable GPU acceleration")
    print(f"   → Command: pip3 install pycuda")

print("\n" + "=" * 80)
print("Test Complete!")
print("=" * 80)
