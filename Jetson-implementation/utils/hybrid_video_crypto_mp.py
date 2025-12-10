"""
Multi-Processing Hybrid Video Encryption Module
Uses multiprocessing instead of threading to bypass Python's GIL
True parallel execution across multiple CPU cores
"""

import numpy as np
import hashlib
import sys
import os
from multiprocessing import Pool, cpu_count
import functools

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


# Global worker function for encryption (must be at module level for pickling)
def _encrypt_channel_worker(args):
    """
    Worker function for channel encryption.
    Must be at module level to be picklable by multiprocessing.
    """
    channel_data, perm, keystream, num_rounds = args
    result = channel_data.copy()
    
    for rnd in range(num_rounds):
        # Confusion: permutation
        result = result[perm[:len(result)]]
        
        # Diffusion: feedback encryption with keystream rotation
        offset = rnd * 13
        ks = np.roll(keystream, -offset)[:len(result)]
        result = feedback_encrypt_channel(result, ks)
    
    return result


def _decrypt_channel_worker(args):
    """
    Worker function for channel decryption.
    Must be at module level to be picklable by multiprocessing.
    """
    channel_data, inv_perm, keystream, num_rounds = args
    result = channel_data.copy()
    
    for rnd in range(num_rounds - 1, -1, -1):
        # Reverse diffusion
        offset = rnd * 13
        ks = np.roll(keystream, -offset)[:len(result)]
        result = feedback_decrypt_channel(result, ks)
        
        # Reverse confusion: inverse permutation
        result = result[inv_perm[:len(result)]]
    
    return result


class HybridVideoEncryptionMP:
    """
    Multi-processing hybrid video encryption using 4 chaotic maps.
    Uses multiprocessing.Pool to bypass Python's GIL.
    True parallel execution across multiple CPU cores.
    """
    
    def __init__(self, frame_width=640, frame_height=480,
                 secret_key="default_secret_key_2025",
                 # Lorenz parameters
                 lorenz_sigma=10.0, lorenz_beta=8/3, lorenz_rho=28.0,
                 # Rossler parameters
                 rossler_a=0.2, rossler_b=0.2, rossler_c=5.7,
                 # Henon parameters
                 henon_a=1.4, henon_b=0.3,
                 # Tent parameters
                 tent_mu=1.9999,
                 # Integration settings
                 dt=0.01, t_end=8.0,
                 # Hash block size
                 hash_block_size=128,
                 # Number of rounds
                 num_rounds=1,
                 # Processing parameters
                 num_processes=3):
        """
        Initialize multi-processing hybrid video encryption.
        
        Args:
            secret_key: Secret key for encryption (string)
            num_processes: Number of worker processes (3 for RGB parallelism)
        """
        self.width = frame_width
        self.height = frame_height
        self.secret_key = secret_key
        
        # Derive initial conditions from secret key
        initial_conditions = self._derive_initial_conditions(secret_key)
        
        # Store parameters with derived initial conditions
        self.lorenz_params = (lorenz_sigma, lorenz_beta, lorenz_rho, 
                             initial_conditions['lorenz_x0'], 
                             initial_conditions['lorenz_y0'], 
                             initial_conditions['lorenz_z0'], 
                             dt, t_end)
        self.rossler_params = (rossler_a, rossler_b, rossler_c,
                              initial_conditions['rossler_x0'], 
                              initial_conditions['rossler_y0'], 
                              initial_conditions['rossler_z0'], 
                              dt, t_end)
        self.henon_params = (henon_a, henon_b, 
                            initial_conditions['henon_x0'], 
                            initial_conditions['henon_y0'])
        self.tent_params = (initial_conditions['tent_x0'], tent_mu)
        
        self.hash_block_size = hash_block_size
        self.num_rounds = num_rounds
        self.num_processes = num_processes
        
        # Create process pool
        self.pool = Pool(processes=num_processes)
        
        # Pre-generate sequences
        self.n_pixels = frame_width * frame_height
        self._generate_all_sequences()
        self._prepare_keystreams()
    
    def _derive_initial_conditions(self, secret_key):
        """
        Derive chaotic map initial conditions from secret key using SHA-256.
        Same implementation as HybridVideoEncryption for consistency.
        """
        key_bytes = secret_key.encode('utf-8')
        
        hashes = []
        for i in range(10):
            h = hashlib.sha256(key_bytes + str(i).encode()).digest()
            value = int.from_bytes(h[:8], byteorder='big') / (2**64)
            hashes.append(value)
        
        lorenz_x0 = hashes[0] * 2 - 1
        lorenz_y0 = hashes[1] * 2 - 1
        lorenz_z0 = hashes[2] * 2 - 1
        rossler_x0 = hashes[3] * 2 - 1
        rossler_y0 = hashes[4] * 2 - 1
        rossler_z0 = hashes[5] * 2 - 1
        henon_x0 = hashes[6] * 2 - 1
        henon_y0 = hashes[7] * 2 - 1
        tent_x0 = 0.1 + hashes[8] * 0.8
        
        return {
            'lorenz_x0': lorenz_x0, 'lorenz_y0': lorenz_y0, 'lorenz_z0': lorenz_z0,
            'rossler_x0': rossler_x0, 'rossler_y0': rossler_y0, 'rossler_z0': rossler_z0,
            'henon_x0': henon_x0, 'henon_y0': henon_y0, 'tent_x0': tent_x0
        }
        
    def _generate_all_sequences(self):
        """Generate all chaotic sequences."""
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
        
    def encrypt_frame(self, frame):
        """
        Encrypt single video frame using multi-processing.
        RGB channels are processed in parallel across worker processes.
        
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
        
        # Prepare arguments for worker processes
        args_list = [
            (r_vec, self.perm_r, self.ks_r, self.num_rounds),
            (g_vec, self.perm_g, self.ks_g, self.num_rounds),
            (b_vec, self.perm_b, self.ks_b, self.num_rounds)
        ]
        
        # Process channels in parallel
        results = self.pool.map(_encrypt_channel_worker, args_list)
        
        # Reconstruct frame
        encrypted = np.zeros((H, W, C), dtype=np.uint8)
        encrypted[:, :, 0] = results[0].reshape(H, W)
        encrypted[:, :, 1] = results[1].reshape(H, W)
        encrypted[:, :, 2] = results[2].reshape(H, W)
        
        return encrypted
        
    def decrypt_frame(self, encrypted_frame):
        """
        Decrypt single video frame using multi-processing.
        RGB channels are processed in parallel across worker processes.
        
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
        
        # Prepare arguments for worker processes
        args_list = [
            (c_r, self.inv_perm_r, self.ks_r, self.num_rounds),
            (c_g, self.inv_perm_g, self.ks_g, self.num_rounds),
            (c_b, self.inv_perm_b, self.ks_b, self.num_rounds)
        ]
        
        # Process channels in parallel
        results = self.pool.map(_decrypt_channel_worker, args_list)
        
        # Reconstruct frame
        decrypted = np.zeros((H, W, C), dtype=np.uint8)
        decrypted[:, :, 0] = results[0].reshape(H, W)
        decrypted[:, :, 1] = results[1].reshape(H, W)
        decrypted[:, :, 2] = results[2].reshape(H, W)
        
        return decrypted
        
    def update_key(self, seed):
        """
        Update encryption key based on seed.
        
        Args:
            seed: seed value for key regeneration
        """
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
        """Cleanup process pool on destruction."""
        if hasattr(self, 'pool'):
            self.pool.close()
            self.pool.join()
