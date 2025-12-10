#!/usr/bin/env python3
"""
Offline Video Encryptor
Encrypts a video file and saves encrypted version to disk
No real-time streaming - just file-to-file encryption
"""

import cv2
import numpy as np
import time
import argparse
import pickle
import os
import sys

# Add parent directory to path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

from utils.hybrid_video_crypto import HybridVideoEncryption
from utils.hybrid_video_crypto_mp import HybridVideoEncryptionMP
from utils.video_crypto import FastVideoEncryption, LightweightVideoEncryption


def encrypt_video_file(input_path, output_path, mode='hybrid', width=320, height=240, use_threads=False):
    """
    Encrypt video file and save to disk.
    
    Args:
        input_path: Input video file path
        output_path: Output encrypted file path
        mode: Encryption mode (hybrid, fast, lightweight)
        width: Frame width
        height: Frame height
        use_threads: Use multi-threaded encryption (for hybrid mode)
    """
    print("\n" + "="*60)
    print("OFFLINE VIDEO ENCRYPTION")
    print("="*60)
    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")
    print(f"Mode:   {mode}{'_MT' if use_threads and mode == 'hybrid' else ''}")
    print(f"Size:   {width}×{height}")
    if use_threads and mode == 'hybrid':
        print(f"Threads: Multi-threaded (3 workers for RGB)")
    
    # Open input video
    cap = cv2.VideoCapture(input_path)
    if not cap.isOpened():
        print(f"Error: Could not open {input_path}")
        return False
    
    # Get video properties
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    fps = cap.get(cv2.CAP_PROP_FPS)
    
    print(f"Frames: {total_frames}")
    print(f"FPS:    {fps:.2f}")
    
    # Initialize encryptor
    print("\nInitializing encryptor...")
    init_start = time.time()
    
    if mode == 'hybrid':
        if use_threads:
            encryptor = HybridVideoEncryptionMP(width, height, num_processes=3)
            print("Using multi-processing (3 processes) for true parallelism")
        else:
            encryptor = HybridVideoEncryption(width, height)
    elif mode == 'fast':
        encryptor = FastVideoEncryption(width, height)
    else:
        encryptor = LightweightVideoEncryption(width, height)
    
    init_time = time.time() - init_start
    print(f"✓ Initialization: {init_time:.2f}s")
    
    # Prepare output
    encrypted_data = {
        'metadata': {
            'mode': mode,
            'width': width,
            'height': height,
            'fps': fps,
            'total_frames': total_frames,
            'original_file': os.path.basename(input_path)
        },
        'frames': []
    }
    
    # Encrypt all frames
    print(f"\nEncrypting {total_frames} frames...")
    start_time = time.time()
    frame_count = 0
    enc_times = []
    
    while True:
        ret, frame = cap.read()
        if not ret:
            break
        
        # Resize if needed
        if frame.shape[1] != width or frame.shape[0] != height:
            frame = cv2.resize(frame, (width, height))
        
        # Convert BGR to RGB
        frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        
        # Encrypt
        enc_start = time.time()
        encrypted_frame = encryptor.encrypt_frame(frame_rgb)
        enc_time = time.time() - enc_start
        enc_times.append(enc_time)
        
        # Save encrypted frame
        encrypted_data['frames'].append(encrypted_frame)
        
        frame_count += 1
        
        # Progress
        if frame_count % 30 == 0 or frame_count == total_frames:
            elapsed = time.time() - start_time
            fps_current = frame_count / elapsed if elapsed > 0 else 0
            avg_enc = np.mean(enc_times) * 1000
            enc_fps = 1000 / avg_enc if avg_enc > 0 else 0
            progress = (frame_count / total_frames) * 100
            print(f"  [{progress:5.1f}%] Frame {frame_count}/{total_frames} | "
                  f"Enc: {avg_enc:.2f}ms ({enc_fps:.1f} FPS) | Overall: {fps_current:.1f} FPS")
    
    cap.release()
    
    # Save to file
    print(f"\nSaving encrypted data to: {output_path}")
    save_start = time.time()
    
    with open(output_path, 'wb') as f:
        pickle.dump(encrypted_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    save_time = time.time() - save_start
    file_size = os.path.getsize(output_path) / (1024 * 1024)  # MB
    
    # Summary
    total_time = time.time() - start_time
    avg_fps = frame_count / total_time if total_time > 0 else 0
    
    print(f"\n" + "="*60)
    print("ENCRYPTION COMPLETE")
    print("="*60)
    print(f"Total frames:      {frame_count}")
    print(f"Total time:        {total_time:.2f}s")
    print(f"Average FPS:       {avg_fps:.2f}")
    print(f"Avg encrypt time:  {np.mean(enc_times)*1000:.2f}ms")
    print(f"File size:         {file_size:.2f} MB")
    print(f"Save time:         {save_time:.2f}s")
    print(f"Output file:       {output_path}")
    print("="*60)
    
    return True


def main():
    parser = argparse.ArgumentParser(description='Offline Video Encryption')
    parser.add_argument('--input', '-i', required=True,
                       help='Input video file')
    parser.add_argument('--output', '-o', required=True,
                       help='Output encrypted file (.enc)')
    parser.add_argument('--mode', '-m', choices=['hybrid', 'fast', 'lightweight'],
                       default='hybrid',
                       help='Encryption mode (default: hybrid)')
    parser.add_argument('--width', '-w', type=int, default=320,
                       help='Frame width (default: 320)')
    parser.add_argument('--height', '-H', type=int, default=240,
                       help='Frame height (default: 240)')
    parser.add_argument('--threads', '-t', action='store_true',
                       help='Use multi-processing for 2-3x speedup (hybrid mode only)')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        return
    
    # Encrypt video
    success = encrypt_video_file(
        args.input,
        args.output,
        mode=args.mode,
        width=args.width,
        height=args.height,
        use_threads=args.threads
    )
    
    if success:
        print("\n✓ Video encrypted successfully!")
        print(f"\nTo decrypt and play:")
        print(f"  python3 decrypt_video_file.py --input {args.output} --mode {args.mode}")


if __name__ == '__main__':
    main()
