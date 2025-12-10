"""
Shared Keystream Manager for TX/RX on Same Machine
Avoids duplicate sequence generation when TX and RX run together
"""

import os
import pickle
import numpy as np

KEYSTREAM_CACHE_DIR = "/tmp/video_encryption_cache"

class SharedKeystreamManager:
    """
    Manages shared keystreams between TX and RX to avoid duplicate work.
    When TX and RX run on same machine, RX can reuse TX's sequences.
    """
    
    def __init__(self, width=640, height=480, mode='hybrid'):
        self.width = width
        self.height = height
        self.mode = mode
        self.cache_file = f"{KEYSTREAM_CACHE_DIR}/keystream_{mode}_{width}x{height}.pkl"
        
    def save_keystreams(self, encryptor):
        """
        Save TX keystreams to shared cache.
        Called by TX after initialization.
        
        Args:
            encryptor: Initialized encryption object
        """
        os.makedirs(KEYSTREAM_CACHE_DIR, exist_ok=True)
        
        cache_data = {
            'width': self.width,
            'height': self.height,
            'mode': self.mode,
        }
        
        # Save based on mode
        if self.mode == 'hybrid':
            cache_data.update({
                'lorenz_seq': encryptor.lorenz_seq,
                'rossler_seq': encryptor.rossler_seq,
                'henon_seq': encryptor.henon_seq,
                'tent_seq': encryptor.tent_seq,
                'ks_r': encryptor.ks_r,
                'ks_g': encryptor.ks_g,
                'ks_b': encryptor.ks_b,
                'perm_r': encryptor.perm_r,
                'perm_g': encryptor.perm_g,
                'perm_b': encryptor.perm_b,
                'inv_perm_r': encryptor.inv_perm_r,
                'inv_perm_g': encryptor.inv_perm_g,
                'inv_perm_b': encryptor.inv_perm_b,
            })
        elif self.mode == 'fast':
            cache_data.update({
                'keystream': encryptor.keystream,
                'perm_indices': encryptor.perm_indices,
                'inv_perm_indices': encryptor.inv_perm_indices,
            })
        else:  # lightweight
            cache_data.update({
                'keystream': encryptor.keystream,
            })
        
        with open(self.cache_file, 'wb') as f:
            pickle.dump(cache_data, f)
        
        print(f"✓ Keystreams cached to: {self.cache_file}")
        
    def load_keystreams(self, decryptor):
        """
        Load cached keystreams into RX decryptor.
        Called by RX before starting reception.
        
        Args:
            decryptor: Decryptor object to populate
            
        Returns:
            True if loaded, False if cache not found
        """
        if not os.path.exists(self.cache_file):
            return False
            
        try:
            with open(self.cache_file, 'rb') as f:
                cache_data = pickle.load(f)
            
            # Verify dimensions match
            if (cache_data['width'] != self.width or 
                cache_data['height'] != self.height or
                cache_data['mode'] != self.mode):
                print("⚠ Cache mismatch, regenerating...")
                return False
            
            # Load based on mode
            if self.mode == 'hybrid':
                decryptor.lorenz_seq = cache_data['lorenz_seq']
                decryptor.rossler_seq = cache_data['rossler_seq']
                decryptor.henon_seq = cache_data['henon_seq']
                decryptor.tent_seq = cache_data['tent_seq']
                decryptor.ks_r = cache_data['ks_r']
                decryptor.ks_g = cache_data['ks_g']
                decryptor.ks_b = cache_data['ks_b']
                decryptor.perm_r = cache_data['perm_r']
                decryptor.perm_g = cache_data['perm_g']
                decryptor.perm_b = cache_data['perm_b']
                decryptor.inv_perm_r = cache_data['inv_perm_r']
                decryptor.inv_perm_g = cache_data['inv_perm_g']
                decryptor.inv_perm_b = cache_data['inv_perm_b']
            elif self.mode == 'fast':
                decryptor.keystream = cache_data['keystream']
                decryptor.perm_indices = cache_data['perm_indices']
                decryptor.inv_perm_indices = cache_data['inv_perm_indices']
            else:  # lightweight
                decryptor.keystream = cache_data['keystream']
            
            print(f"✓ Keystreams loaded from cache (skipped generation)")
            return True
            
        except Exception as e:
            print(f"⚠ Cache load failed: {e}")
            return False
    
    def clear_cache(self):
        """Remove cached keystreams."""
        if os.path.exists(self.cache_file):
            os.remove(self.cache_file)
            print(f"✓ Cache cleared: {self.cache_file}")

def clear_all_cache():
    """Clear all cached keystreams."""
    if os.path.exists(KEYSTREAM_CACHE_DIR):
        import shutil
        shutil.rmtree(KEYSTREAM_CACHE_DIR)
        print(f"✓ All caches cleared")
