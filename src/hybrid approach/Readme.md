# Hybrid Chaos-Based Image Encryption (Improved Method)

This folder contains the MATLAB implementation of an **improved hybrid chaos-based image encryption algorithm** that combines multiple chaotic systems to enhance security, diffusion strength, and statistical randomness.

The main implementation file is:

- **`hybridapproachmain.m`**

---

## Overview of the Hybrid Method

This hybrid encryption approach integrates **continuous-time and discrete-time chaotic systems** to construct a strong encryption framework with multiple layers of security.

### Chaotic Systems Used
The algorithm combines the following chaotic maps:

- **Lorenz system** – used as a continuous-time entropy source
- **Rossler system** – provides additional nonlinear chaotic dynamics
- **Henon map** – discrete chaotic map for permutation and mixing
- **Tent map** – discrete chaotic map for fast randomness generation

Each chaotic system contributes independently to the key generation, permutation, and diffusion processes.

---

## Key Features of the Implementation

The proposed hybrid method includes the following improvements:

- **Per-channel keystream generation** (independent R, G, B processing)
- **SHA-256–based whitening** of chaotic sequences for enhanced randomness
- **Per-channel pixel permutation (confusion stage)**
- **Feedback-based diffusion (CBC-like structure)**
- **Two rounds of confusion + diffusion**
- **Key sensitivity due to multiple chaotic parameters**
- **Exact decryption with correct key reproduction**

---

## Encryption Process (High-Level)

1. Generate chaotic sequences using:
   - Lorenz system (ODE-based)
   - Rossler system (ODE-based)
   - Henon map
   - Tent map
2. Normalize and mix chaotic outputs to create per-channel chaotic signals.
3. Apply **SHA-256 hashing** to whiten chaotic data into byte-level keystreams.
4. Generate **channel-wise permutation indices**.
5. Perform **two rounds** of:
   - Pixel permutation (confusion)
   - Feedback diffusion using XOR-like operations
6. Reconstruct encrypted RGB image.

---

## Decryption Process

- Uses the **same chaotic parameters and initial conditions**.
- Performs inverse operations in reverse order:
  - Reverse diffusion
  - Inverse permutation
- Perfect reconstruction is achieved if keys are correct.

---

## Security and Evaluation Metrics

The script evaluates encryption strength using:

- **Histogram analysis**
- **Chi-square statistical test**
- **Pixel correlation analysis (Horizontal, Vertical, Diagonal)**
- **Encryption time measurement**

These metrics demonstrate resistance against:
- Statistical attacks
- Correlation-based attacks
- Differential key sensitivity attacks



