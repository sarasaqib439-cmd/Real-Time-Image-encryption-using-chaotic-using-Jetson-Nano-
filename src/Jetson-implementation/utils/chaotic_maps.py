"""
Chaotic Map Generators for Image Encryption
Optimized for Jetson Nano with NumPy vectorization
"""

import numpy as np
from scipy.integrate import odeint


def generate_lorenz_map(sigma=10.0, beta=8/3, rho=28.0, 
                        x0=0.1, y0=0.0, z0=0.1, 
                        dt=0.01, t_end=200.0):
    """
    Generate Lorenz chaotic attractor sequence.
    
    Args:
        sigma, beta, rho: Lorenz system parameters
        x0, y0, z0: Initial conditions
        dt: Time step
        t_end: End time for integration
        
    Returns:
        numpy array of shape (n_points, 3) containing x, y, z trajectories
    """
    def lorenz_system(state, t):
        x, y, z = state
        dx = sigma * (y - x)
        dy = x * (rho - z) - y
        dz = x * y - beta * z
        return [dx, dy, dz]
    
    t = np.arange(0, t_end, dt)
    initial_state = [x0, y0, z0]
    trajectory = odeint(lorenz_system, initial_state, t)
    
    return trajectory


def generate_rossler_map(a=0.2, b=0.2, c=5.7,
                         x0=0.1, y0=0.1, z0=0.1,
                         dt=0.01, t_end=200.0):
    """
    Generate Rossler chaotic attractor sequence.
    
    Args:
        a, b, c: Rossler system parameters
        x0, y0, z0: Initial conditions
        dt: Time step
        t_end: End time for integration
        
    Returns:
        numpy array of shape (n_points, 3) containing x, y, z trajectories
    """
    def rossler_system(state, t):
        x, y, z = state
        dx = -y - z
        dy = x + a * y
        dz = b + z * (x - c)
        return [dx, dy, dz]
    
    t = np.arange(0, t_end, dt)
    initial_state = [x0, y0, z0]
    trajectory = odeint(rossler_system, initial_state, t)
    
    return trajectory


def generate_henon_map(a=1.4, b=0.3, x0=0.1, y0=0.1, n_iterations=100000):
    """
    Generate Henon map chaotic sequence.
    
    Args:
        a, b: Henon map parameters
        x0, y0: Initial conditions
        n_iterations: Number of iterations
        
    Returns:
        numpy array of shape (n_iterations, 2) containing x, y sequences
    """
    sequence = np.zeros((n_iterations, 2))
    x, y = x0, y0
    
    for i in range(n_iterations):
        x_next = 1 - a * x**2 + y
        y_next = b * x
        sequence[i] = [x_next, y_next]
        x, y = x_next, y_next
    
    return sequence


def generate_tent_map(mu=1.9999, x0=0.12345, n_iterations=100000):
    """
    Generate Tent map chaotic sequence.
    
    Args:
        mu: Tent map parameter (typically close to 2)
        x0: Initial condition (0 < x0 < 1)
        n_iterations: Number of iterations
        
    Returns:
        numpy array of shape (n_iterations,) containing the sequence
    """
    sequence = np.zeros(n_iterations)
    x = x0
    
    for i in range(n_iterations):
        if x < 0.5:
            x = mu * x
        else:
            x = mu * (1 - x)
        sequence[i] = x
    
    return sequence


def normalize_sequence(seq):
    """
    Normalize sequence to [0, 1] range.
    
    Args:
        seq: Input sequence (numpy array)
        
    Returns:
        Normalized sequence in [0, 1] range
    """
    min_val = np.min(seq)
    max_val = np.max(seq)
    return (seq - min_val) / (max_val - min_val + 1e-10)


def extend_sequence(seq, required_length):
    """
    Extend sequence to required length by repeating if necessary.
    
    Args:
        seq: Input sequence
        required_length: Desired length
        
    Returns:
        Extended sequence
    """
    if len(seq) >= required_length:
        return seq[:required_length]
    
    repetitions = (required_length // len(seq)) + 1
    extended = np.tile(seq, repetitions)
    return extended[:required_length]
