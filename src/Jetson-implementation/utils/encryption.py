"""
Encryption and Decryption Functions
Includes confusion, diffusion, and feedback mechanisms
Optimized for Jetson Nano
"""

import numpy as np
import hashlib


def sha256_hash_doubles(data):
    """
    Generate SHA-256 hash from double array.
    
    Args:
        data: numpy array of floats
        
    Returns:
        32 bytes as uint8 numpy array
    """
    # Convert doubles to bytes
    byte_data = data.tobytes()
    # Compute SHA-256 hash
    hash_obj = hashlib.sha256(byte_data)
    digest = hash_obj.digest()
    # Convert to uint8 numpy array
    return np.frombuffer(digest, dtype=np.uint8)


def generate_whitened_keystream(source_vec, n_needed, block_size=32):
    """
    Generate whitened keystream using SHA-256 hashing.
    
    Args:
        source_vec: Normalized chaotic sequence
        n_needed: Number of bytes needed
        block_size: Number of doubles per hash block
        
    Returns:
        uint8 keystream array of length n_needed
    """
    keystream = np.zeros(n_needed, dtype=np.uint8)
    pos = 0
    idx = 0
    L = len(source_vec)
    
    while pos < n_needed:
        # Extract block with wraparound
        start_idx = (idx * block_size) % L
        block_indices = np.arange(start_idx, start_idx + block_size) % L
        block = source_vec[block_indices].astype(np.float64)
        
        # Hash the block
        hash_bytes = sha256_hash_doubles(block)
        
        # Copy to keystream
        take = min(len(hash_bytes), n_needed - pos)
        keystream[pos:pos + take] = hash_bytes[:take]
        pos += take
        idx += 1
    
    return keystream


def feedback_encrypt_channel(plaintext, keystream):
    """
    Encrypt single channel with feedback (CBC-like).
    
    Args:
        plaintext: uint8 array
        keystream: uint8 array (same length)
        
    Returns:
        Encrypted uint8 array
    """
    n = len(plaintext)
    ciphertext = np.zeros(n, dtype=np.uint8)
    
    # Use first keystream byte as IV
    iv = keystream[0]
    
    # First element
    ciphertext[0] = plaintext[0] ^ keystream[0] ^ iv
    
    # Remaining elements with feedback
    for i in range(1, n):
        ciphertext[i] = plaintext[i] ^ keystream[i] ^ ciphertext[i-1]
    
    return ciphertext


def feedback_decrypt_channel(ciphertext, keystream):
    """
    Decrypt single channel with feedback (CBC-like).
    
    Args:
        ciphertext: uint8 array
        keystream: uint8 array (same length)
        
    Returns:
        Decrypted uint8 array
    """
    n = len(ciphertext)
    plaintext = np.zeros(n, dtype=np.uint8)
    
    # Use first keystream byte as IV
    iv = keystream[0]
    
    # First element
    plaintext[0] = ciphertext[0] ^ keystream[0] ^ iv
    
    # Remaining elements with feedback
    for i in range(1, n):
        plaintext[i] = ciphertext[i] ^ keystream[i] ^ ciphertext[i-1]
    
    return plaintext


def confusion_stage(image, chaotic_seq):
    """
    Confusion stage - permutation based on chaotic sequence.
    
    Args:
        image: Input image (H, W, C) uint8
        chaotic_seq: Chaotic sequence for generating permutation
        
    Returns:
        Confused image
    """
    H, W, C = image.shape
    N = H * W * C
    
    # Flatten image
    flat_img = image.flatten()
    
    # Generate permutation indices from chaotic sequence
    # Normalize and convert to indices
    normalized = (chaotic_seq[:N] - np.min(chaotic_seq[:N])) / \
                 (np.max(chaotic_seq[:N]) - np.min(chaotic_seq[:N]) + 1e-10)
    
    # Use argsort to create permutation
    perm_indices = np.argsort(normalized)
    
    # Apply permutation
    confused_flat = flat_img[perm_indices]
    
    # Reshape back
    confused_img = confused_flat.reshape(H, W, C)
    
    return confused_img.astype(np.uint8)


def reverse_confusion(confused_img, chaotic_seq):
    """
    Reverse confusion stage - inverse permutation.
    
    Args:
        confused_img: Confused image (H, W, C) uint8
        chaotic_seq: Same chaotic sequence used for encryption
        
    Returns:
        Original image
    """
    H, W, C = confused_img.shape
    N = H * W * C
    
    # Flatten image
    flat_img = confused_img.flatten()
    
    # Generate same permutation indices
    normalized = (chaotic_seq[:N] - np.min(chaotic_seq[:N])) / \
                 (np.max(chaotic_seq[:N]) - np.min(chaotic_seq[:N]) + 1e-10)
    perm_indices = np.argsort(normalized)
    
    # Create inverse permutation
    inv_perm = np.argsort(perm_indices)
    
    # Apply inverse permutation
    original_flat = flat_img[inv_perm]
    
    # Reshape back
    original_img = original_flat.reshape(H, W, C)
    
    return original_img.astype(np.uint8)


def diffusion_stage(confused_img, chaotic_seq):
    """
    Diffusion stage - XOR with chaotic sequence.
    
    Args:
        confused_img: Confused image (H, W, C) uint8
        chaotic_seq: Chaotic sequence for XOR operation
        
    Returns:
        Encrypted image
    """
    H, W, C = confused_img.shape
    N = H * W * C
    
    # Flatten image
    flat_img = confused_img.flatten()
    
    # Generate keystream from chaotic sequence
    normalized = (chaotic_seq[:N] - np.min(chaotic_seq[:N])) / \
                 (np.max(chaotic_seq[:N]) - np.min(chaotic_seq[:N]) + 1e-10)
    keystream = (normalized * 255).astype(np.uint8)
    
    # XOR diffusion
    encrypted_flat = flat_img ^ keystream
    
    # Reshape back
    encrypted_img = encrypted_flat.reshape(H, W, C)
    
    return encrypted_img.astype(np.uint8)


def reverse_diffusion(encrypted_img, chaotic_seq):
    """
    Reverse diffusion stage - XOR again (symmetric operation).
    
    Args:
        encrypted_img: Encrypted image (H, W, C) uint8
        chaotic_seq: Same chaotic sequence used for encryption
        
    Returns:
        Confused image (before encryption)
    """
    # XOR is its own inverse
    return diffusion_stage(encrypted_img, chaotic_seq)


def permute_channel(channel_data, perm_indices):
    """
    Permute a single channel using given indices.
    
    Args:
        channel_data: 1D array of channel data
        perm_indices: Permutation indices
        
    Returns:
        Permuted channel data
    """
    return channel_data[perm_indices]


def inverse_permute_channel(permuted_data, perm_indices):
    """
    Inverse permutation of a single channel.
    
    Args:
        permuted_data: 1D array of permuted data
        perm_indices: Original permutation indices
        
    Returns:
        Original channel data
    """
    inv_perm = np.argsort(perm_indices)
    return permuted_data[inv_perm]
