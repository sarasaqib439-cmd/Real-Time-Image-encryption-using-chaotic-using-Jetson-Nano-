"""
Socket Communication Utilities for Video Streaming
Optimized for low-latency transmission
"""

import socket
import struct
import pickle
import numpy as np
import cv2


class VideoSocket:
    """Handle socket communication for video frames."""
    
    def __init__(self, host='localhost', port=9999):
        self.host = host
        self.port = port
        self.sock = None
        self.conn = None
        self.buffer_size = 32768  # 32KB buffer (Jetson optimized: reduced from 64KB)
        
    def create_server(self):
        """Create server socket (for TX)."""
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.bind((self.host, self.port))
        self.sock.listen(5)
        print(f"Server listening on {self.host}:{self.port}")
        
    def accept_connection(self):
        """Accept incoming connection."""
        self.conn, addr = self.sock.accept()
        print(f"Connected to {addr}")
        return self.conn
        
    def connect_to_server(self):
        """Connect to server (for RX)."""
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((self.host, self.port))
        print(f"Connected to server at {self.host}:{self.port}")
        
    def send_frame(self, frame, metadata=None):
        """
        Send frame over socket with size header.
        
        Args:
            frame: numpy array (image)
            metadata: optional dictionary with frame info
        """
        # Prepare data
        data = {
            'frame': frame,
            'shape': frame.shape,
            'dtype': str(frame.dtype)
        }
        if metadata:
            data['metadata'] = metadata
            
        # Serialize
        serialized = pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL)
        size = len(serialized)
        
        # Send size header (4 bytes) then data
        try:
            if self.conn:
                self.conn.sendall(struct.pack('>I', size))
                self.conn.sendall(serialized)
            else:
                self.sock.sendall(struct.pack('>I', size))
                self.sock.sendall(serialized)
            return True
        except (BrokenPipeError, ConnectionResetError):
            print("Connection lost")
            return False
            
    def receive_frame(self):
        """
        Receive frame from socket.
        
        Returns:
            frame, metadata or None, None if error
        """
        try:
            # Receive size header (4 bytes)
            size_data = self._recv_all(4)
            if not size_data:
                return None, None
                
            size = struct.unpack('>I', size_data)[0]
            
            # Receive frame data
            data = self._recv_all(size)
            if not data:
                return None, None
                
            # Deserialize
            unpacked = pickle.loads(data)
            frame = unpacked['frame']
            metadata = unpacked.get('metadata', {})
            
            return frame, metadata
            
        except Exception as e:
            print(f"Error receiving frame: {e}")
            return None, None
            
    def _recv_all(self, n):
        """Helper to receive n bytes."""
        data = bytearray()
        while len(data) < n:
            packet = self.sock.recv(n - len(data))
            if not packet:
                return None
            data.extend(packet)
        return bytes(data)
        
    def close(self):
        """Close socket connections."""
        if self.conn:
            self.conn.close()
        if self.sock:
            self.sock.close()
        print("Socket closed")


class FrameBuffer:
    """Simple frame buffer for smooth playback."""
    
    def __init__(self, max_size=10):
        self.buffer = []
        self.max_size = max_size
        
    def add(self, frame, metadata=None):
        """Add frame to buffer."""
        if len(self.buffer) >= self.max_size:
            self.buffer.pop(0)
        self.buffer.append((frame, metadata))
        
    def get(self):
        """Get frame from buffer."""
        if self.buffer:
            return self.buffer.pop(0)
        return None, None
        
    def is_empty(self):
        """Check if buffer is empty."""
        return len(self.buffer) == 0
        
    def size(self):
        """Get current buffer size."""
        return len(self.buffer)


def compress_frame(frame, quality=80):
    """
    Compress frame using JPEG encoding.
    
    Args:
        frame: numpy array (image)
        quality: JPEG quality (0-100)
        
    Returns:
        compressed bytes
    """
    encode_param = [int(cv2.IMWRITE_JPEG_QUALITY), quality]
    result, encoded = cv2.imencode('.jpg', frame, encode_param)
    return encoded.tobytes() if result else None


def decompress_frame(data):
    """
    Decompress JPEG frame.
    
    Args:
        data: compressed bytes
        
    Returns:
        numpy array (image)
    """
    nparr = np.frombuffer(data, np.uint8)
    frame = cv2.imdecode(nparr, cv2.IMREAD_COLOR)
    return frame


def calculate_fps(frame_times, window=30):
    """
    Calculate FPS from frame timestamps.
    
    Args:
        frame_times: list of timestamps
        window: number of frames to average
        
    Returns:
        current FPS
    """
    if len(frame_times) < 2:
        return 0.0
        
    recent = frame_times[-min(window, len(frame_times)):]
    if len(recent) < 2:
        return 0.0
        
    time_diff = recent[-1] - recent[0]
    if time_diff > 0:
        return (len(recent) - 1) / time_diff
    return 0.0


def resize_frame(frame, target_width=640):
    """
    Resize frame maintaining aspect ratio.
    
    Args:
        frame: input frame
        target_width: desired width
        
    Returns:
        resized frame
    """
    h, w = frame.shape[:2]
    aspect_ratio = w / h
    target_height = int(target_width / aspect_ratio)
    return cv2.resize(frame, (target_width, target_height))
