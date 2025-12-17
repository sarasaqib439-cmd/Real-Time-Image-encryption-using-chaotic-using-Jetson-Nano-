"""
Multi-Threaded Hybrid Video Encryption Module
Uses all 4 chaotic maps (Lorenz, Rossler, Henon, Tent) with parallel channel processing
Optimized for Jetson Nano's 4-core ARM processor
"""

import numpy as np
import hashlib
import sys
import os
from concurrent.futures import ThreadPoolExecutor
import threading

# Add parent directories to path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '../..'))
sys.path.insert(0, parent_dir)

try:
    from utils.chaotic_maps import (
        generate_lorenz_map,
        generate_rossler_map,
        generate_henon_map,
        generate_tent_map,
        normalize_sequence
    )
    from utils.encryption import (
        generate_whitened_keystream,
        feedback_encrypt_channel,
        feedback_decrypt_channel
    )
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(current_dir, '../../..')))
    from python_implementation.utils.chaotic_maps import (
        generate_lorenz_map,
        generate_rossler_map,
        generate_henon_map,
        generate_tent_map,
        normalize_sequence
    )
    from python_implementation.utils.encryption import (
        generate_whitened_keystream,
        feedback_encrypt_channel,
        feedback_decrypt_channel
    )


class HybridVideoEncryptionMT:
    """
    Multi-threaded hybrid video encryption using 4 chaotic maps.
    Lorenz + Rossler + Henon + Tent with SHA-256 whitening.
    Parallelizes RGB channel processing across multiple CPU cores.
    """
    
    def __init__(self, frame_width=640, frame_height=480,
                 # Lorenz parameters
                 lorenz_sigma=10.0, lorenz_beta=8/3, lorenz_rho=28.0,
                 lorenz_x0=0.1, lorenz_y0=0.0, lorenz_z0=0.1,
                 # Rossler parameters
                 rossler_a=0.2, rossler_b=0.2, rossler_c=5.7,
                 rossler_x0=0.1, rossler_y0=0.1, rossler_z0=0.1,
                 # Henon parameters
                 henon_a=1.4, henon_b=0.3,
                 henon_x0=0.1, henon_y0=0.1,
                 # Tent parameters
                 tent_x0=0.12345, tent_mu=1.9999,
                 # Integration settings (Jetson Nano optimized)
                 dt=0.01, t_end=8.0,
                 # Hash block size (Jetson Nano optimized)
                 hash_block_size=128,
                 # Number of rounds
                 num_rounds=1,
                 # Threading parameters
                 max_workers=3):  # 3 workers for RGB channels
        """
        Initialize multi-threaded hybrid video encryption.
        
        Args:
            frame_width, frame_height: Video dimensions
            lorenz_*, rossler_*, henon_*, tent_*: Chaotic map parameters
            dt, t_end: Integration parameters for continuous systems
            hash_block_size: SHA-256 block size for whitening
            num_rounds: Encryption rounds (1 for real-time)
            max_workers: Number of worker threads (3 for RGB parallelism)
        """
        self.width = frame_width
        self.height = frame_height
        
        # Store parameters
        self.lorenz_params = (lorenz_sigma, lorenz_beta, lorenz_rho, 
                             lorenz_x0, lorenz_y0, lorenz_z0, dt, t_end)
        self.rossler_params = (rossler_a, rossler_b, rossler_c,
                              rossler_x0, rossler_y0, rossler_z0, dt, t_end)
        self.henon_params = (henon_a, henon_b, henon_x0, henon_y0)
        self.tent_params = (tent_x0, tent_mu)
        
        self.hash_block_size = hash_block_size
        self.num_rounds = num_rounds
        self.max_workers = max_workers
        
        # Thread pool for parallel processing
        self.executor = ThreadPoolExecutor(max_workers=max_workers)
        self.lock = threading.Lock()
        
        # Pre-generate sequences
        self.n_pixels = frame_width * frame_height
        self._generate_all_sequences()
        self._prepare_keystreams()
        
    def _generate_all_sequences(self):
        """Generate all chaotic sequences (Jetson Nano optimized)."""
        length_needed = self.n_pixels + 1000
        
        # Lorenz
        self.lorenz_seq = generate_lorenz_map(*self.lorenz_params)
        
        # Rossler
        self.rossler_seq = generate_rossler_map(*self.rossler_params)
        
        # Henon
        self.henon_seq = generate_henon_map(*self.henon_params, length_needed)
        
        # Tent
        self.tent_seq = generate_tent_map(*self.tent_params, length_needed)
        
        # Extend Lorenz/Rossler if needed
        if len(self.lorenz_seq) < length_needed:
            rep_factor = (length_needed // len(self.lorenz_seq)) + 2
            self.lorenz_seq = np.tile(self.lorenz_seq, (rep_factor, 1))
        self.lorenz_seq = self.lorenz_seq[:length_needed]
        
        if len(self.rossler_seq) < length_needed:
            rep_factor = (length_needed // len(self.rossler_seq)) + 2
            self.rossler_seq = np.tile(self.rossler_seq, (rep_factor, 1))
        self.rossler_seq = self.rossler_seq[:length_needed]
        
    def _prepare_keystreams(self):
        """Prepare per-channel keystreams and permutations."""
        n_chan = self.n_pixels
        
        # Normalize components
        a1 = normalize_sequence(self.lorenz_seq[:, 0])
        a2 = normalize_sequence(self.rossler_seq[:, 0])
        a3 = normalize_sequence(self.henon_seq[:, 0])
        a4 = normalize_sequence(self.tent_seq)
        
        # Create channel-specific mixes
        mix_r = normalize_sequence(0.6 * a1 + 0.4 * a3)
        mix_g = normalize_sequence(0.5 * a2 + 0.5 * a4)
        mix_b = normalize_sequence(0.5 * (a1 + a2))
        
        # Generate whitened keystreams
        self.ks_r = generate_whitened_keystream(mix_r, n_chan, self.hash_block_size)
        self.ks_g = generate_whitened_keystream(mix_g, n_chan, self.hash_block_size)
        self.ks_b = generate_whitened_keystream(mix_b, n_chan, self.hash_block_size)
        
        # Create per-channel permutations
        temp_r = (a1[:n_chan] * 1e12).astype(np.uint64)
        temp_g = (a2[:n_chan] * 1e12).astype(np.uint64)
        temp_b = (a3[:n_chan] * 1e12).astype(np.uint64)
        
        perm_seed_r = ((temp_r + (temp_g << 5)) % (2**52)).astype(float)
        perm_seed_g = ((temp_g + (temp_b << 7)) % (2**52)).astype(float)
        perm_seed_b = ((temp_b + (temp_r << 11)) % (2**52)).astype(float)
        
        self.perm_r = np.argsort(perm_seed_r)
        self.perm_g = np.argsort(perm_seed_g)
        self.perm_b = np.argsort(perm_seed_b)
        
        # Inverse permutations for decryption
        self.inv_perm_r = np.argsort(self.perm_r)
        self.inv_perm_g = np.argsort(self.perm_g)
        self.inv_perm_b = np.argsort(self.perm_b)
    
    def _encrypt_channel(self, channel_data, perm, keystream):
        """
        Encrypt a single channel (called by worker thread).
        
        Args:
            channel_data: 1D numpy array of channel pixels
            perm: Permutation array for this channel
            keystream: Keystream for this channel
            
        Returns:
            Encrypted channel data
        """
        result = channel_data.copy()
        
        # Multi-round encryption
        for rnd in range(self.num_rounds):
            # Confusion: permutation
            result = result[perm[:len(result)]]
            
            # Diffusion: feedback encryption with keystream rotation
            offset = rnd * 13
            ks = np.roll(keystream, -offset)[:len(result)]
            result = feedback_encrypt_channel(result, ks)
        
        return result
    
    def _decrypt_channel(self, channel_data, inv_perm, keystream):
        """
        Decrypt a single channel (called by worker thread).
        
        Args:
            channel_data: 1D numpy array of encrypted channel pixels
            inv_perm: Inverse permutation array for this channel
            keystream: Keystream for this channel
            
        Returns:
            Decrypted channel data
        """
        result = channel_data.copy()
        
        # Multi-round decryption (reverse order)
        for rnd in range(self.num_rounds - 1, -1, -1):
            # Reverse diffusion
            offset = rnd * 13
            ks = np.roll(keystream, -offset)[:len(result)]
            result = feedback_decrypt_channel(result, ks)
            
            # Reverse confusion: inverse permutation
            result = result[inv_perm[:len(result)]]
        
        return result
        
    def encrypt_frame(self, frame):
        """
        Encrypt single video frame using multi-threaded hybrid approach.
        RGB channels are processed in parallel across 3 worker threads.
        
        Args:
            frame: numpy array (H, W, 3) uint8
            
        Returns:
            encrypted frame
        """
        H, W, C = frame.shape
        
        # Flatten channels
        r_vec = frame[:, :, 0].flatten()
        g_vec = frame[:, :, 1].flatten()
        b_vec = frame[:, :, 2].flatten()
        
        # Submit channel encryption tasks to thread pool
        future_r = self.executor.submit(self._encrypt_channel, r_vec, self.perm_r, self.ks_r)
        future_g = self.executor.submit(self._encrypt_channel, g_vec, self.perm_g, self.ks_g)
        future_b = self.executor.submit(self._encrypt_channel, b_vec, self.perm_b, self.ks_b)
        
        # Wait for all channels to complete
        p_r = future_r.result()
        p_g = future_g.result()
        p_b = future_b.result()
        
        # Reconstruct frame
        encrypted = np.zeros((H, W, C), dtype=np.uint8)
        encrypted[:, :, 0] = p_r.reshape(H, W)
        encrypted[:, :, 1] = p_g.reshape(H, W)
        encrypted[:, :, 2] = p_b.reshape(H, W)
        
        return encrypted
        
    def decrypt_frame(self, encrypted_frame):
        """
        Decrypt single video frame using multi-threaded approach.
        RGB channels are processed in parallel across 3 worker threads.
        
        Args:
            encrypted_frame: encrypted numpy array (H, W, 3) uint8
            
        Returns:
            decrypted frame
        """
        H, W, C = encrypted_frame.shape
        
        # Flatten channels
        c_r = encrypted_frame[:, :, 0].flatten()
        c_g = encrypted_frame[:, :, 1].flatten()
        c_b = encrypted_frame[:, :, 2].flatten()
        
        # Submit channel decryption tasks to thread pool
        future_r = self.executor.submit(self._decrypt_channel, c_r, self.inv_perm_r, self.ks_r)
        future_g = self.executor.submit(self._decrypt_channel, c_g, self.inv_perm_g, self.ks_g)
        future_b = self.executor.submit(self._decrypt_channel, c_b, self.inv_perm_b, self.ks_b)
        
        # Wait for all channels to complete
        d_r = future_r.result()
        d_g = future_g.result()
        d_b = future_b.result()
        
        # Reconstruct frame
        decrypted = np.zeros((H, W, C), dtype=np.uint8)
        decrypted[:, :, 0] = d_r.reshape(H, W)
        decrypted[:, :, 1] = d_g.reshape(H, W)
        decrypted[:, :, 2] = d_b.reshape(H, W)
        
        return decrypted
        
    def update_key(self, seed):
        """
        Update encryption key based on seed.
        Allows key to change per frame if needed.
        
        Args:
            seed: seed value for key regeneration
        """
        with self.lock:  # Thread-safe key update
            # Hash seed to get new initial conditions
            hash_val = hashlib.sha256(str(seed).encode()).digest()
            
            # Update all initial conditions
            self.lorenz_params = (
                self.lorenz_params[0], self.lorenz_params[1], self.lorenz_params[2],
                (int.from_bytes(hash_val[:4], 'big') % 10000) / 10000,
                (int.from_bytes(hash_val[4:8], 'big') % 10000) / 10000,
                (int.from_bytes(hash_val[8:12], 'big') % 10000) / 10000,
                self.lorenz_params[6], self.lorenz_params[7]
            )
            
            self.rossler_params = (
                self.rossler_params[0], self.rossler_params[1], self.rossler_params[2],
                (int.from_bytes(hash_val[12:16], 'big') % 10000) / 10000,
                (int.from_bytes(hash_val[16:20], 'big') % 10000) / 10000,
                (int.from_bytes(hash_val[20:24], 'big') % 10000) / 10000,
                self.rossler_params[6], self.rossler_params[7]
            )
            
            # Regenerate sequences and keystreams
            self._generate_all_sequences()
            self._prepare_keystreams()
    
    def __del__(self):
        """Cleanup thread pool on destruction."""
        if hasattr(self, 'executor'):
            self.executor.shutdown(wait=True)
