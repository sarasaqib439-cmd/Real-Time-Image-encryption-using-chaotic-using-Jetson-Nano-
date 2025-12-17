"""
Analysis Tools for Encrypted Images
Including correlation analysis, chi-square test, and histogram plotting
Optimized for Jetson Nano
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr


def chi_square_histogram(image, method_name='Image'):
    """
    Calculate chi-square values for each channel's histogram.
    
    Args:
        image: Input RGB image (H, W, 3) uint8
        method_name: Name of the encryption method
        
    Returns:
        Array of chi-square values [R, G, B]
    """
    chi2_values = np.zeros(3)
    
    for ch in range(3):
        channel_data = image[:, :, ch].flatten()
        N = len(channel_data)
        
        # Calculate histogram
        counts, _ = np.histogram(channel_data, bins=256, range=(0, 256))
        
        # Expected count per bin (uniform distribution)
        E = N / 256.0
        
        # Chi-square formula: sum((observed - expected)^2 / expected)
        chi2_values[ch] = np.sum((counts - E)**2 / E)
    
    # Print results
    print(f'\nChi-Square Values ({method_name}):')
    print(f'Red Channel   : {chi2_values[0]:.4f}')
    print(f'Green Channel : {chi2_values[1]:.4f}')
    print(f'Blue Channel  : {chi2_values[2]:.4f}')
    
    return chi2_values


def plot_histogram(image, title='Image Histogram', save_path=None):
    """
    Plot histogram for all three channels.
    
    Args:
        image: Input RGB image (H, W, 3) uint8
        title: Title for the plot
        save_path: Optional path to save the figure
    """
    colors = ['red', 'green', 'blue']
    channel_names = ['Red', 'Green', 'Blue']
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    for ch in range(3):
        channel_data = image[:, :, ch].flatten()
        axes[ch].hist(channel_data, bins=256, range=(0, 256), 
                     color=colors[ch], alpha=0.7, edgecolor='black')
        axes[ch].set_title(f'{channel_names[ch]} Channel')
        axes[ch].set_xlabel('Pixel Value')
        axes[ch].set_ylabel('Frequency')
        axes[ch].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Histogram saved to {save_path}')
    else:
        plt.show()
    
    plt.close()


def correlation_analysis(image, title='Correlation Analysis', save_path=None, sample_size=5000):
    """
    Perform correlation analysis on adjacent pixels in three directions.
    
    Args:
        image: Input RGB image (H, W, 3) uint8
        title: Title for the plot
        save_path: Optional path to save the figure
        sample_size: Number of pixels to sample for scatter plots
    """
    channel_names = ['Red', 'Green', 'Blue']
    colors = ['red', 'green', 'blue']
    directions = ['Horizontal', 'Vertical', 'Diagonal']
    
    # Store correlation coefficients
    results = np.zeros((3, 3))  # 3 channels x 3 directions
    
    fig, axes = plt.subplots(3, 3, figsize=(14, 12))
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    for ch in range(3):
        channel = image[:, :, ch].astype(float)
        
        # Horizontal correlation
        x_h = channel[:, :-1].flatten()
        y_h = channel[:, 1:].flatten()
        if len(x_h) > sample_size:
            indices = np.random.choice(len(x_h), sample_size, replace=False)
            x_h_sample, y_h_sample = x_h[indices], y_h[indices]
        else:
            x_h_sample, y_h_sample = x_h, y_h
        
        axes[ch, 0].scatter(x_h_sample, y_h_sample, s=1, c=colors[ch], alpha=0.5)
        axes[ch, 0].set_title(f'{channel_names[ch]} - Horizontal')
        axes[ch, 0].set_xlabel('Pixel(i,j)')
        axes[ch, 0].set_ylabel('Pixel(i,j+1)')
        axes[ch, 0].grid(True, alpha=0.3)
        
        corr_h, _ = pearsonr(x_h, y_h)
        results[ch, 0] = corr_h
        
        # Vertical correlation
        x_v = channel[:-1, :].flatten()
        y_v = channel[1:, :].flatten()
        if len(x_v) > sample_size:
            indices = np.random.choice(len(x_v), sample_size, replace=False)
            x_v_sample, y_v_sample = x_v[indices], y_v[indices]
        else:
            x_v_sample, y_v_sample = x_v, y_v
        
        axes[ch, 1].scatter(x_v_sample, y_v_sample, s=1, c=colors[ch], alpha=0.5)
        axes[ch, 1].set_title(f'{channel_names[ch]} - Vertical')
        axes[ch, 1].set_xlabel('Pixel(i,j)')
        axes[ch, 1].set_ylabel('Pixel(i+1,j)')
        axes[ch, 1].grid(True, alpha=0.3)
        
        corr_v, _ = pearsonr(x_v, y_v)
        results[ch, 1] = corr_v
        
        # Diagonal correlation
        x_d = channel[:-1, :-1].flatten()
        y_d = channel[1:, 1:].flatten()
        if len(x_d) > sample_size:
            indices = np.random.choice(len(x_d), sample_size, replace=False)
            x_d_sample, y_d_sample = x_d[indices], y_d[indices]
        else:
            x_d_sample, y_d_sample = x_d, y_d
        
        axes[ch, 2].scatter(x_d_sample, y_d_sample, s=1, c=colors[ch], alpha=0.5)
        axes[ch, 2].set_title(f'{channel_names[ch]} - Diagonal')
        axes[ch, 2].set_xlabel('Pixel(i,j)')
        axes[ch, 2].set_ylabel('Pixel(i+1,j+1)')
        axes[ch, 2].grid(True, alpha=0.3)
        
        corr_d, _ = pearsonr(x_d, y_d)
        results[ch, 2] = corr_d
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Correlation plot saved to {save_path}')
    else:
        plt.show()
    
    plt.close()
    
    # Print correlation summary
    print(f'\n================= Correlation Summary: {title} =================')
    print('              Horizontal     Vertical       Diagonal')
    print('---------------------------------------------------------------')
    print(f'Red     :     {results[0, 0]:0.4f}        {results[0, 1]:0.4f}         {results[0, 2]:0.4f}')
    print(f'Green   :     {results[1, 0]:0.4f}        {results[1, 1]:0.4f}         {results[1, 2]:0.4f}')
    print(f'Blue    :     {results[2, 0]:0.4f}        {results[2, 1]:0.4f}         {results[2, 2]:0.4f}')
    print('===============================================================\n')
    
    return results


def display_images(original, encrypted, decrypted=None, save_path=None):
    """
    Display original, encrypted, and optionally decrypted images.
    
    Args:
        original: Original image
        encrypted: Encrypted image
        decrypted: Decrypted image (optional)
        save_path: Optional path to save the figure
    """
    if decrypted is not None:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        axes[0].imshow(original)
        axes[0].set_title('Original Image')
        axes[0].axis('off')
        
        axes[1].imshow(encrypted)
        axes[1].set_title('Encrypted Image')
        axes[1].axis('off')
        
        axes[2].imshow(decrypted)
        axes[2].set_title('Decrypted Image')
        axes[2].axis('off')
    else:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        axes[0].imshow(original)
        axes[0].set_title('Original Image')
        axes[0].axis('off')
        
        axes[1].imshow(encrypted)
        axes[1].set_title('Encrypted Image')
        axes[1].axis('off')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Image comparison saved to {save_path}')
    else:
        plt.show()
    
    plt.close()


def calculate_entropy(image):
    """
    Calculate information entropy for each channel.
    
    Args:
        image: Input RGB image (H, W, 3) uint8
        
    Returns:
        Array of entropy values [R, G, B]
    """
    entropy_values = np.zeros(3)
    
    for ch in range(3):
        channel_data = image[:, :, ch].flatten()
        
        # Calculate histogram
        hist, _ = np.histogram(channel_data, bins=256, range=(0, 256))
        
        # Calculate probabilities (remove zeros to avoid log(0))
        hist = hist[hist > 0]
        probs = hist / np.sum(hist)
        
        # Calculate entropy: -sum(p * log2(p))
        entropy_values[ch] = -np.sum(probs * np.log2(probs))
    
    print(f'\nEntropy Values:')
    print(f'Red Channel   : {entropy_values[0]:.6f}')
    print(f'Green Channel : {entropy_values[1]:.6f}')
    print(f'Blue Channel  : {entropy_values[2]:.6f}')
    
    return entropy_values


def calculate_npcr_uaci(image1, image2):
    """
    Calculate NPCR (Number of Pixels Change Rate) and 
    UACI (Unified Average Changing Intensity) between two images.
    
    Args:
        image1: First image (H, W, 3) uint8
        image2: Second image (H, W, 3) uint8
        
    Returns:
        Tuple of (NPCR, UACI) for each channel
    """
    H, W, C = image1.shape
    total_pixels = H * W
    
    npcr_values = np.zeros(3)
    uaci_values = np.zeros(3)
    
    for ch in range(3):
        ch1 = image1[:, :, ch].astype(float)
        ch2 = image2[:, :, ch].astype(float)
        
        # NPCR: percentage of different pixels
        diff_pixels = np.sum(ch1 != ch2)
        npcr_values[ch] = (diff_pixels / total_pixels) * 100
        
        # UACI: average intensity difference
        uaci_values[ch] = (np.sum(np.abs(ch1 - ch2)) / (total_pixels * 255)) * 100
    
    print(f'\nNPCR Values (%):')
    print(f'Red Channel   : {npcr_values[0]:.4f}')
    print(f'Green Channel : {npcr_values[1]:.4f}')
    print(f'Blue Channel  : {npcr_values[2]:.4f}')
    
    print(f'\nUACI Values (%):')
    print(f'Red Channel   : {uaci_values[0]:.4f}')
    print(f'Green Channel : {uaci_values[1]:.4f}')
    print(f'Blue Channel  : {uaci_values[2]:.4f}')
    
    return npcr_values, uaci_values
