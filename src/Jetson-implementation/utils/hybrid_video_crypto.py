"""
Hybrid Video Encryption Module
Uses all 4 chaotic maps (Lorenz, Rossler, Henon, Tent) for stronger security
"""

import numpy as np
import hashlib
import sys
import os

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


class HybridVideoEncryption:
    """
    Hybrid video encryption using 4 chaotic maps.
    Lorenz + Rossler + Henon + Tent with SHA-256 whitening.
    Optimized for real-time performance.
    Supports secret key for deriving initial conditions.
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
                 # Integration settings (Jetson Nano optimized)
                 dt=0.01, t_end=8.0,  # Reduced to 8.0 for Jetson Nano (30% less memory)
                 # Hash block size (Jetson Nano optimized)
                 hash_block_size=128,  # Larger block = 2× fewer hash ops (faster on Jetson)
                 # Number of rounds
                 num_rounds=1):  # Single round for real-time on Jetson Nano
        """
        Initialize hybrid video encryption.
        
        Args:
            frame_width, frame_height: Video dimensions
            secret_key: Secret key for encryption (string)
            lorenz_*, rossler_*, henon_*, tent_*: Chaotic map parameters
            dt, t_end: Integration parameters for continuous systems
            hash_block_size: SHA-256 block size for whitening
            num_rounds: Encryption rounds (1 for real-time, 2+ for security)
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
        
        # Pre-generate sequences
        self.n_pixels = frame_width * frame_height
        self._generate_all_sequences()
        self._prepare_keystreams()
    
    def _derive_initial_conditions(self, secret_key):
        """
        Derive chaotic map initial conditions from secret key using SHA-256.
        This makes the encryption key-dependent and secure.
        
        Args:
            secret_key: User's secret key (string)
            
        Returns:
            Dictionary of initial conditions for all chaotic maps
        """
        # Hash the key multiple times to get different values
        key_bytes = secret_key.encode('utf-8')
        
        # Generate 10 different hashes for 10 initial values
        hashes = []
        for i in range(10):
            h = hashlib.sha256(key_bytes + str(i).encode()).digest()
            # Convert first 8 bytes to float in range (0, 1)
            value = int.from_bytes(h[:8], byteorder='big') / (2**64)
            # Scale to appropriate range for each chaotic map
            hashes.append(value)
        
        # Lorenz: x, y, z in range [-1, 1]
        lorenz_x0 = hashes[0] * 2 - 1
        lorenz_y0 = hashes[1] * 2 - 1
        lorenz_z0 = hashes[2] * 2 - 1
        
        # Rossler: x, y, z in range [-1, 1]
        rossler_x0 = hashes[3] * 2 - 1
        rossler_y0 = hashes[4] * 2 - 1
        rossler_z0 = hashes[5] * 2 - 1
        
        # Henon: x, y in range [-1, 1]
        henon_x0 = hashes[6] * 2 - 1
        henon_y0 = hashes[7] * 2 - 1
        
        # Tent: x in range (0, 1)
        tent_x0 = 0.1 + hashes[8] * 0.8  # Avoid extremes
        
        return {
            'lorenz_x0': lorenz_x0,
            'lorenz_y0': lorenz_y0,
            'lorenz_z0': lorenz_z0,
            'rossler_x0': rossler_x0,
            'rossler_y0': rossler_y0,
            'rossler_z0': rossler_z0,
            'henon_x0': henon_x0,
            'henon_y0': henon_y0,
            'tent_x0': tent_x0
        }
        
    def _generate_all_sequences(self):
        """Generate all chaotic sequences (Jetson Nano optimized)."""
        # Reduced buffer from 2000 to 1000 for Jetson memory savings
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
        Encrypt single video frame using hybrid approach.
        
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
        
        p_r, p_g, p_b = r_vec, g_vec, b_vec
        
        # Multi-round encryption
        for rnd in range(self.num_rounds):
            # Confusion: per-channel permutation
            p_r = p_r[self.perm_r[:len(p_r)]]
            p_g = p_g[self.perm_g[:len(p_g)]]
            p_b = p_b[self.perm_b[:len(p_b)]]
            
            # Diffusion: feedback encryption with keystream rotation
            offset = rnd * 13
            kr = np.roll(self.ks_r, -offset)[:len(p_r)]
            kg = np.roll(self.ks_g, -offset)[:len(p_g)]
            kb = np.roll(self.ks_b, -offset)[:len(p_b)]
            
            p_r = feedback_encrypt_channel(p_r, kr)
            p_g = feedback_encrypt_channel(p_g, kg)
            p_b = feedback_encrypt_channel(p_b, kb)
        
        # Reconstruct frame
        encrypted = np.zeros((H, W, C), dtype=np.uint8)
        encrypted[:, :, 0] = p_r.reshape(H, W)
        encrypted[:, :, 1] = p_g.reshape(H, W)
        encrypted[:, :, 2] = p_b.reshape(H, W)
        
        return encrypted
        
    def decrypt_frame(self, encrypted_frame):
        """
        Decrypt single video frame.
        
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
        
        # Multi-round decryption (reverse order)
        for rnd in range(self.num_rounds - 1, -1, -1):
            # Reverse diffusion
            offset = rnd * 13
            kr = np.roll(self.ks_r, -offset)[:len(c_r)]
            kg = np.roll(self.ks_g, -offset)[:len(c_g)]
            kb = np.roll(self.ks_b, -offset)[:len(c_b)]
            
            c_r = feedback_decrypt_channel(c_r, kr)
            c_g = feedback_decrypt_channel(c_g, kg)
            c_b = feedback_decrypt_channel(c_b, kb)
            
            # Reverse confusion: inverse permutation
            c_r = c_r[self.inv_perm_r[:len(c_r)]]
            c_g = c_g[self.inv_perm_g[:len(c_g)]]
            c_b = c_b[self.inv_perm_b[:len(c_b)]]
        
        # Reconstruct frame
        decrypted = np.zeros((H, W, C), dtype=np.uint8)
        decrypted[:, :, 0] = c_r.reshape(H, W)
        decrypted[:, :, 1] = c_g.reshape(H, W)
        decrypted[:, :, 2] = c_b.reshape(H, W)
        
        return decrypted
        
    def update_key(self, seed):
        """
        Update encryption key based on seed.
        Allows key to change per frame if needed.
        
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
