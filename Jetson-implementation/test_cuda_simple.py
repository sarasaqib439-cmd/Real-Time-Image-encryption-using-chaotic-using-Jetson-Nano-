#!/usr/bin/env python3
"""
Simple CUDA Test - Verify PyCUDA works and test basic performance
Run this first to make sure CUDA is working before full benchmark
"""

import numpy as np
import time

print("=" * 70)
print("CUDA Quick Test - PyCUDA Verification")
print("=" * 70)

# Step 1: Check if PyCUDA is available
print("\n[Step 1/4] Checking PyCUDA installation...")
try:
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    print("✅ PyCUDA imported successfully")
except ImportError as e:
    print(f"❌ PyCUDA import failed: {e}")
    print("\nInstall with: pip3 install pycuda")
    exit(1)

# Step 2: Get CUDA device info
print("\n[Step 2/4] Getting CUDA device information...")
try:
    device = cuda.Device(0)
    print(f"✅ CUDA Device: {device.name()}")
    print(f"   Compute Capability: {device.compute_capability()}")
    print(f"   Total Memory: {device.total_memory() / 1024**2:.0f} MB")
    print(f"   Multiprocessors: {device.get_attribute(cuda.device_attribute.MULTIPROCESSOR_COUNT)}")
except Exception as e:
    print(f"❌ Failed to get device info: {e}")
    exit(1)

# Step 3: Test simple CUDA kernel
print("\n[Step 3/4] Testing simple CUDA kernel compilation...")

simple_kernel = """
__global__ void test_xor_kernel(unsigned char *input, unsigned char *output, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        output[idx] = input[idx] ^ 0xFF;
    }
}
"""

try:
    mod = SourceModule(simple_kernel)
    test_func = mod.get_function("test_xor_kernel")
    print("✅ CUDA kernel compiled successfully")
except Exception as e:
    print(f"❌ Kernel compilation failed: {e}")
    exit(1)

# Step 4: Test kernel execution
print("\n[Step 4/4] Testing CUDA kernel execution...")

try:
    # Create test data
    n = 1024
    test_data = np.random.randint(0, 256, n, dtype=np.uint8)
    output = np.zeros(n, dtype=np.uint8)
    
    # Allocate GPU memory
    d_input = cuda.mem_alloc(test_data.nbytes)
    d_output = cuda.mem_alloc(output.nbytes)
    
    # Copy to GPU
    cuda.memcpy_htod(d_input, test_data)
    
    # Run kernel
    block_size = 256
    grid_size = (n + block_size - 1) // block_size
    test_func(d_input, d_output, np.int32(n),
              block=(block_size, 1, 1),
              grid=(grid_size, 1))
    
    # Copy back
    cuda.memcpy_dtoh(output, d_output)
    
    # Verify
    expected = test_data ^ 0xFF
    if np.array_equal(output, expected):
        print("✅ CUDA kernel execution successful - results verified!")
    else:
        print("❌ CUDA kernel produced incorrect results")
        exit(1)
        
except Exception as e:
    print(f"❌ Kernel execution failed: {e}")
    exit(1)

# Step 5: Simple performance test
print("\n" + "=" * 70)
print("Performance Test - XOR Operation")
print("=" * 70)

# Test different sizes
sizes = [
    (76800, "320×240"),    # 320x240 pixels
    (307200, "640×480"),   # 640x480 pixels
]

for size, name in sizes:
    print(f"\n{name} ({size} bytes)")
    
    test_data = np.random.randint(0, 256, size, dtype=np.uint8)
    
    # CPU version
    start = time.time()
    for _ in range(100):
        cpu_result = test_data ^ 0xFF
    cpu_time = time.time() - start
    cpu_rate = (size * 100) / cpu_time / 1e6  # MB/s
    
    # GPU version
    d_input = cuda.mem_alloc(test_data.nbytes)
    d_output = cuda.mem_alloc(test_data.nbytes)
    output = np.empty_like(test_data)
    
    block_size = 256
    grid_size = (size + block_size - 1) // block_size
    
    # Warm up
    for _ in range(10):
        cuda.memcpy_htod(d_input, test_data)
        test_func(d_input, d_output, np.int32(size),
                  block=(block_size, 1, 1), grid=(grid_size, 1))
        cuda.memcpy_dtoh(output, d_output)
    
    # Benchmark
    start = time.time()
    for _ in range(100):
        cuda.memcpy_htod(d_input, test_data)
        test_func(d_input, d_output, np.int32(size),
                  block=(block_size, 1, 1), grid=(grid_size, 1))
        cuda.memcpy_dtoh(output, d_output)
    gpu_time = time.time() - start
    gpu_rate = (size * 100) / gpu_time / 1e6  # MB/s
    
    speedup = cpu_time / gpu_time
    
    print(f"   CPU: {cpu_time*10:.2f} ms (100 ops), {cpu_rate:.1f} MB/s")
    print(f"   GPU: {gpu_time*10:.2f} ms (100 ops), {gpu_rate:.1f} MB/s")
    print(f"   Speedup: {speedup:.2f}×")

print("\n" + "=" * 70)
print("CUDA Test Complete!")
print("=" * 70)
print("\n✅ PyCUDA is working correctly on your Jetson Nano!")
print("   You can now run the full encryption benchmark.")
print("\nNext steps:")
print("   1. Copy the CUDA implementation files to utils/")
print("   2. Run: python3 benchmark_performance.py")
print("=" * 70)
