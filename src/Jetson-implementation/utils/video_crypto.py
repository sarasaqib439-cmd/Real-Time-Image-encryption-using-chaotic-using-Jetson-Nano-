"""
Real-time Video Encryption Module
Fast chaotic encryption for video frames
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
    from utils.chaotic_maps import generate_henon_map, normalize_sequence
except ImportError:
    # If running from video_encryption directory, try relative import
    sys.path.insert(0, os.path.abspath(os.path.join(current_dir, '../../..')))
    from python_implementation.utils.chaotic_maps import generate_henon_map, normalize_sequence


class FastVideoEncryption:
    """
    Fast video encryption using Henon map.
    Optimized for real-time performance on Jetson Nano.
    """
    
    def __init__(self, frame_width=640, frame_height=480,
                 henon_a=1.4, henon_b=0.3, 
                 henon_x0=0.9, henon_y0=0.5):
        """
        Initialize video encryption.
        
        Args:
            frame_width: Video frame width
            frame_height: Video frame height
            henon_a, henon_b: Henon map parameters
            henon_x0, henon_y0: Initial conditions
        """
        self.width = frame_width
        self.height = frame_height
        self.henon_a = henon_a
        self.henon_b = henon_b
        self.henon_x0 = henon_x0
        self.henon_y0 = henon_y0
        
        # Pre-generate chaotic sequence for one frame
        self.n_pixels = frame_width * frame_height * 3
        self._generate_keystream()
        
        # Pre-compute permutation indices
        self._generate_permutation()
        
    def _generate_keystream(self):
        """Generate chaotic keystream."""
        # Generate Henon sequence
        henon_seq = generate_henon_map(
            self.henon_a, self.henon_b,
            self.henon_x0, self.henon_y0,
            self.n_pixels
        )
        
        # Normalize and convert to keystream
        normalized = normalize_sequence(henon_seq[:, 0])
        self.keystream = (normalized * 255).astype(np.uint8)
        
    def _generate_permutation(self):
        """Generate permutation indices."""
        # Use second column of Henon for permutation
        henon_seq = generate_henon_map(
            self.henon_a, self.henon_b,
            self.henon_x0, self.henon_y0,
            self.n_pixels
        )
        normalized = normalize_sequence(henon_seq[:, 1])
        self.perm_indices = np.argsort(normalized)
        self.inv_perm_indices = np.argsort(self.perm_indices)
        
    def encrypt_frame(self, frame):
        """
        Encrypt single video frame.
        
        Args:
            frame: numpy array (H, W, 3) uint8
            
        Returns:
            encrypted frame
        """
        H, W, C = frame.shape
        
        # Flatten frame
        flat = frame.flatten()
        
        # Confusion: permute pixels
        permuted = flat[self.perm_indices[:len(flat)]]
        
        # Diffusion: XOR with keystream
        encrypted = permuted ^ self.keystream[:len(flat)]
        
        # Reshape back
        return encrypted.reshape(H, W, C).astype(np.uint8)
        
    def decrypt_frame(self, encrypted_frame):
        """
        Decrypt single video frame.
        
        Args:
            encrypted_frame: encrypted numpy array (H, W, 3) uint8
            
        Returns:
            decrypted frame
        """
        H, W, C = encrypted_frame.shape
        
        # Flatten frame
        flat = encrypted_frame.flatten()
        
        # Reverse diffusion: XOR with keystream
        diffused = flat ^ self.keystream[:len(flat)]
        
        # Reverse confusion: inverse permutation
        decrypted = diffused[self.inv_perm_indices[:len(flat)]]
        
        # Reshape back
        return decrypted.reshape(H, W, C).astype(np.uint8)
        
    def update_key(self, seed):
        """
        Update encryption key based on seed.
        Allows key to change per frame if needed.
        
        Args:
            seed: seed value for key regeneration
        """
        # Hash seed to get new initial conditions
        hash_val = hashlib.sha256(str(seed).encode()).digest()
        self.henon_x0 = (int.from_bytes(hash_val[:4], 'big') % 10000) / 10000
        self.henon_y0 = (int.from_bytes(hash_val[4:8], 'big') % 10000) / 10000
        
        # Regenerate keystream and permutation
        self._generate_keystream()
        self._generate_permutation()


class LightweightVideoEncryption:
    """
    Ultra-lightweight video encryption.
    XOR-only for maximum speed (60+ FPS target).
    """
    
    def __init__(self, frame_width=640, frame_height=480, seed=12345):
        """
        Initialize lightweight encryption.
        
        Args:
            frame_width: Video frame width
            frame_height: Video frame height
            seed: Random seed for keystream
        """
        self.width = frame_width
        self.height = frame_height
        self.seed = seed
        
        # Generate simple keystream
        n_pixels = frame_width * frame_height * 3
        np.random.seed(seed)
        self.keystream = np.random.randint(0, 256, n_pixels, dtype=np.uint8)
        
    def encrypt_frame(self, frame):
        """Fast XOR encryption."""
        flat = frame.flatten()
        encrypted = flat ^ self.keystream[:len(flat)]
        return encrypted.reshape(frame.shape).astype(np.uint8)
        
    def decrypt_frame(self, encrypted_frame):
        """Fast XOR decryption (symmetric)."""
        return self.encrypt_frame(encrypted_frame)


class AdaptiveVideoEncryption:
    """
    Adaptive encryption that adjusts to maintain target FPS.
    """
    
    def __init__(self, frame_width=640, frame_height=480, target_fps=30):
        """
        Initialize adaptive encryption.
        
        Args:
            frame_width: Video frame width
            frame_height: Video frame height
            target_fps: Target FPS to maintain
        """
        # Start with fast encryption
        self.fast_enc = FastVideoEncryption(frame_width, frame_height)
        self.light_enc = LightweightVideoEncryption(frame_width, frame_height)
        
        self.target_fps = target_fps
        self.use_fast = True  # Start with full encryption
        self.frame_times = []
        
    def encrypt_frame(self, frame):
        """Encrypt with adaptive algorithm selection."""
        if self.use_fast:
            return self.fast_enc.encrypt_frame(frame)
        else:
            return self.light_enc.encrypt_frame(frame)
            
    def decrypt_frame(self, encrypted_frame):
        """Decrypt with adaptive algorithm selection."""
        if self.use_fast:
            return self.fast_enc.decrypt_frame(encrypted_frame)
        else:
            return self.light_enc.decrypt_frame(encrypted_frame)
            
    def update_performance(self, current_fps):
        """
        Adjust encryption based on performance.
        
        Args:
            current_fps: Current achieved FPS
        """
        # Switch to lightweight if falling below target
        if current_fps < self.target_fps * 0.8:
            self.use_fast = False
            print("Switching to lightweight encryption")
        elif current_fps > self.target_fps * 1.2:
            self.use_fast = True
            print("Switching to full encryption")
