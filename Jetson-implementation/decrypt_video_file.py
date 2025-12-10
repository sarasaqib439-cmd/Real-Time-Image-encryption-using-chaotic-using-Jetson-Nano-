#!/usr/bin/env python3
"""
Offline Video Decryptor
Decrypts an encrypted video file and displays/saves it
No real-time streaming - just file-to-file decryption
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


def decrypt_video_file(input_path, output_path=None, display=True, mode='hybrid', use_threads=False):
    """
    Decrypt video file and display/save.
    
    Args:
        input_path: Encrypted file path
        output_path: Output video file (optional, for saving)
        display: Show video while decrypting
        mode: Encryption mode (must match encryption)
        use_threads: Use multi-threaded decryption (for hybrid mode)
    """
    print("\n" + "="*60)
    print("OFFLINE VIDEO DECRYPTION")
    print("="*60)
    print(f"Input:   {input_path}")
    if output_path:
        print(f"Output:  {output_path}")
    print(f"Display: {'Yes' if display else 'No'}")
    
    # Load encrypted data
    print("\nLoading encrypted data...")
    load_start = time.time()
    
    if not os.path.exists(input_path):
        print(f"Error: File not found: {input_path}")
        return False
    
    with open(input_path, 'rb') as f:
        encrypted_data = pickle.load(f)
    
    load_time = time.time() - load_start
    file_size = os.path.getsize(input_path) / (1024 * 1024)  # MB
    
    # Get metadata
    metadata = encrypted_data['metadata']
    frames = encrypted_data['frames']
    
    print(f"✓ Loaded: {file_size:.2f} MB in {load_time:.2f}s")
    print(f"Mode:     {metadata['mode']}")
    print(f"Size:     {metadata['width']}×{metadata['height']}")
    print(f"FPS:      {metadata['fps']:.2f}")
    print(f"Frames:   {len(frames)}")
    
    # Verify mode matches
    if metadata['mode'] != mode:
        print(f"\n⚠ Warning: Encrypted with '{metadata['mode']}' but decrypting with '{mode}'")
        print(f"  This may produce incorrect results!")
        response = input("Continue anyway? (y/n): ")
        if response.lower() != 'y':
            return False
    
    # Initialize decryptor
    print("\nInitializing decryptor...")
    init_start = time.time()
    
    width = metadata['width']
    height = metadata['height']
    
    if mode == 'hybrid':
        if use_threads:
            decryptor = HybridVideoEncryptionMP(width, height, num_processes=3)
            print(f"Using multi-processing decryption (3 processes)")
        else:
            decryptor = HybridVideoEncryption(width, height)
    elif mode == 'fast':
        decryptor = FastVideoEncryption(width, height)
    else:
        decryptor = LightweightVideoEncryption(width, height)
    
    init_time = time.time() - init_start
    print(f"✓ Initialization: {init_time:.2f}s")
    
    # Setup output video writer if saving
    writer = None
    if output_path:
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        writer = cv2.VideoWriter(output_path, fourcc, metadata['fps'],
                                 (width, height))
        print(f"✓ Output writer initialized")
    
    # Decrypt all frames
    print(f"\nDecrypting {len(frames)} frames...")
    start_time = time.time()
    frame_count = 0
    dec_times = []
    
    for encrypted_frame in frames:
        # Decrypt
        dec_start = time.time()
        decrypted_frame = decryptor.decrypt_frame(encrypted_frame)
        dec_time = time.time() - dec_start
        dec_times.append(dec_time)
        
        # Convert RGB to BGR for display/save
        decrypted_bgr = cv2.cvtColor(decrypted_frame, cv2.COLOR_RGB2BGR)
        
        # Save if output specified
        if writer:
            writer.write(decrypted_bgr)
        
        # Display if requested
        if display:
            cv2.imshow('Decrypted Video', decrypted_bgr)
            
            # Control playback speed
            delay = int(1000 / metadata['fps'])
            if cv2.waitKey(delay) & 0xFF == ord('q'):
                print("\nPlayback stopped by user")
                break
        
        frame_count += 1
        
        # Progress
        if frame_count % 30 == 0 or frame_count == len(frames):
            elapsed = time.time() - start_time
            fps_current = frame_count / elapsed if elapsed > 0 else 0
            avg_dec = np.mean(dec_times) * 1000
            progress = (frame_count / len(frames)) * 100
            print(f"  [{progress:5.1f}%] Frame {frame_count}/{len(frames)} | "
                  f"Dec: {avg_dec:.2f}ms | Speed: {fps_current:.1f} FPS")
    
    # Cleanup
    if writer:
        writer.release()
    if display:
        cv2.destroyAllWindows()
    
    # Summary
    total_time = time.time() - start_time
    avg_fps = frame_count / total_time if total_time > 0 else 0
    
    print(f"\n" + "="*60)
    print("DECRYPTION COMPLETE")
    print("="*60)
    print(f"Total frames:      {frame_count}")
    print(f"Total time:        {total_time:.2f}s")
    print(f"Average FPS:       {avg_fps:.2f}")
    print(f"Avg decrypt time:  {np.mean(dec_times)*1000:.2f}ms")
    if output_path:
        output_size = os.path.getsize(output_path) / (1024 * 1024)
        print(f"Output size:       {output_size:.2f} MB")
        print(f"Output file:       {output_path}")
    print("="*60)
    
    return True


def main():
    parser = argparse.ArgumentParser(description='Offline Video Decryption')
    parser.add_argument('--input', '-i', required=True,
                       help='Input encrypted file (.enc)')
    parser.add_argument('--output', '-o',
                       help='Output video file (optional)')
    parser.add_argument('--mode', '-m', choices=['hybrid', 'fast', 'lightweight'],
                       default='hybrid',
                       help='Decryption mode (must match encryption)')
    parser.add_argument('--no-display', action='store_true',
                       help='Do not display video (only save)')
    parser.add_argument('--threads', '-t', action='store_true',
                       help='Use multi-processing for 2-3x speedup (hybrid mode only)')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        return
    
    # Decrypt video
    success = decrypt_video_file(
        args.input,
        output_path=args.output,
        display=not args.no_display,
        mode=args.mode,
        use_threads=args.threads
    )
    
    if success:
        print("\n✓ Video decrypted successfully!")
        if args.output:
            print(f"\n✓ Saved to: {args.output}")


if __name__ == '__main__':
    main()
