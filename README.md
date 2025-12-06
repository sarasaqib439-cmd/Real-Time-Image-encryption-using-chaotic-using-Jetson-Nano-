# **Chaotic Image Encryption using Multiple Chaotic Maps**

## **Overview**
This project implements and compares chaotic image encryption techniques based on **confusion and diffusion principles**. The focus is to analyze security performance and computational efficiency of different chaotic map–based encryption methods for images.

---

## **Implemented Methods**
- **Henon Map Encryption**
- **Tent Map Encryption**
- **Lorenz–Rossler Chaotic Encryption**
- **Improved Hybrid Chaotic Method (Proposed)**

---

## **Features**
- Pixel-level **confusion and diffusion**
- **Histogram uniformity** after encryption
- **Low pixel correlation** in encrypted images
- **High key sensitivity**
- **Correct decryption** with the same key
- MATLAB-based implementation

---

## **Performance Evaluation**
The following metrics are used for analysis:
- Histogram Analysis  
- Chi-square Test  
- Correlation Analysis (Horizontal, Vertical, Diagonal)  
- Key Sensitivity Test  
- Encryption/Decryption Time  

---

## **Summary of Results**
- Henon map provides fast encryption with moderate security  
- Tent map improves histogram uniformity  
- Lorenz–Rossler method offers stronger security but higher computation time  
- The hybrid method provides a balanced trade-off between security and speed  

---

## **How to Run**
1. Open **MATLAB**
2. Keep the project folder structure unchanged
3. Place the input image in the appropriate directory or update the image path inside the script
4. Run any one of the following main files:

### **Hybrid Method**
- `codes/hybrid approach/hybridapproachmain.m`

### **State-of-the-Art Methods**
- `codes/state of art methods/henonmain.m`
- `codes/state of art methods/tentmain.m`
- `codes/state of art methods/lorenzrosslermain.m`


---

## **Applications**
- Secure image transmission  
- Medical image protection  
- Surveillance data security  
- Research and academic study  
