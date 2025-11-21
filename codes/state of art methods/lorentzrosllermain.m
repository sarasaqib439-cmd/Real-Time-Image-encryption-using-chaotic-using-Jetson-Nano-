%% Chaotic Maps Image Encryption
clc; clear; close all;
%% Read the image
img = imread('C:\Users\rimsha malik\Desktop\PhD NUST\computer vision\assignment 2\project\codes\state of art methods\paper.PNG');
img = uint8(img);

[R, C, ~] = size(img);
N = R*C;

%% Lorenz–Rossler parameters
delta = 20; 
r = 20; 
a = 9; 
beta = 8.5; 
b = 0; 
c = 8;

% initial conditions
x0 = 0.001;  y0 = 0.001;  z0 = 0.1;

%% === 1. Generate chaotic sequences ===
[X, Y, Z] = lorenz_rossler_system(delta, r, a, beta, b, c, [x0; y0; z0], N);

%% === 2. Normalize chaotic sequences ===
Xn = uint8(mod(abs(X)*1e14, 256));
Yn = uint8(mod(abs(Y)*1e14, 256));
Zn = uint8(mod(abs(Z)*1e14, 256));

%% === 3. Split image ===
Rch = img(:,:,1);
Gch = img(:,:,2);
Bch = img(:,:,3);

Rvec = Rch(:);
Gvec = Gch(:);
Bvec = Bch(:);

%% === 4. Sort sequences and obtain confusion index ===
[~, idxX] = sort(Xn);
[~, idxY] = sort(Yn);
[~, idxZ] = sort(Zn);

%% === 5. Confusion (rearranging pixels) ===
shfR = Rvec(idxX);
shfG = Gvec(idxY);
shfB = Bvec(idxZ);

%% === 6. Diffusion (XOR) ===
encR = bitxor(shfR, Xn);
encG = bitxor(shfG, Yn);
encB = bitxor(shfB, Zn);

%% === 7. Merge encrypted channels ===
enc_img = zeros(R,C,3,'uint8');
enc_img(:,:,1) = reshape(encR, [R,C]);
enc_img(:,:,2) = reshape(encG, [R,C]);
enc_img(:,:,3) = reshape(encB, [R,C]);

%% Display
figure;
subplot(1,2,1), imshow(img), title('Original');
subplot(1,2,2), imshow(enc_img), title('Encrypted (Lorenz–Rossler)');

%% %% Histogram Analysis for Lorentz Rossler Map Encryption
fprintf('=== Lorentz Rossler Map Encryption Histogram Analysis ===\n');
% Histogram Analysis
plot_hist(img, 'Original Image');
plot_hist(enc_img, 'Encrypted Image  (Lorenz–Rossler)');
chi_Image   = chi_square_hist(img, 'Original Image');
chi_LR    = chi_square_hist(enc_img,'Lorenz–Rossler' );

%% === DECRYPTION ===
% Input: enc_img (encrypted image)
% Output: dec_img (decrypted image)

[R, C, ~] = size(enc_img);
N = R * C;

%% 1. Regenerate chaotic sequences (must be identical)
[X, Y, Z] = lorenz_rossler_system(delta, r, a, beta, b, c, [x0; y0; z0], N);

%% 2. Normalize sequences (same as encryption)
Xn = uint8(mod(abs(X)*1e14, 256));
Yn = uint8(mod(abs(Y)*1e14, 256));
Zn = uint8(mod(abs(Z)*1e14, 256));

%% 3. Split encrypted image into channels
encR = enc_img(:,:,1); encG = enc_img(:,:,2); encB = enc_img(:,:,3);

encRvec = encR(:); encGvec = encG(:); encBvec = encB(:);

%% 4. Undo diffusion (XOR again with same sequences)
shfR = bitxor(encRvec, Xn);   % reverses encR = bitxor(shfR, Xn)
shfG = bitxor(encGvec, Yn);
shfB = bitxor(encBvec, Zn);

%% 5. Undo confusion (reorder pixels back to original)
[~, idxX] = sort(Xn);   % same index used in encryption
[~, idxY] = sort(Yn);
[~, idxZ] = sort(Zn);

invR = zeros(N,1,'uint8');  invG = zeros(N,1,'uint8');  invB = zeros(N,1,'uint8');

invR(idxX) = shfR;   % place pixels back using the index
invG(idxY) = shfG;
invB(idxZ) = shfB;

%% 6. Merge channels to reconstruct original image
dec_img = zeros(R,C,3,'uint8');
dec_img(:,:,1) = reshape(invR, [R,C]);
dec_img(:,:,2) = reshape(invG, [R,C]);
dec_img(:,:,3) = reshape(invB, [R,C]);

%% 7. Display
figure;
subplot(1,3,1), imshow(enc_img), title('Encrypted Image');
subplot(1,3,2), imshow(dec_img), title('Decrypted Image');
subplot(1,3,3), imshow(img), title('Original Image');

% Optional check
if isequal(uint8(img), dec_img)
    disp('Decryption successful: original image recovered.');
else
    disp('Decryption failed: mismatch detected.');
end


%% Correlation Analysis:
correlation_analysis(img, 'Original Image Correlation');
correlation_analysis(enc_img, '(Lorenz–Rossler) Encrypted Image Correlation');

%% %% Timing Analysis for Lorenz–Rossler Encryption
fprintf('=== Lorenz–Rossler Encryption Timing Analysis ===\n');

% Start timer
tic;

% 1. Generate chaotic sequences
[X, Y, Z] = lorenz_rossler_system(delta, r, a, beta, b, c, [x0; y0; z0], N);

% 2. Normalize chaotic sequences
Xn = uint8(mod(abs(X)*1e14, 256));
Yn = uint8(mod(abs(Y)*1e14, 256));
Zn = uint8(mod(abs(Z)*1e14, 256));

% 3. Split image
Rch = img(:,:,1); Gch = img(:,:,2); Bch = img(:,:,3);
Rvec = Rch(:); Gvec = Gch(:); Bvec = Bch(:);

% 4. Confusion
[~, idxX] = sort(Xn); [~, idxY] = sort(Yn); [~, idxZ] = sort(Zn);
shfR = Rvec(idxX); shfG = Gvec(idxY); shfB = Bvec(idxZ);

% 5. Diffusion (XOR)
encR = bitxor(shfR, Xn);
encG = bitxor(shfG, Yn);
encB = bitxor(shfB, Zn);

% 6. Merge channels
enc_img = zeros(R,C,3,'uint8');
enc_img(:,:,1) = reshape(encR, [R,C]);
enc_img(:,:,2) = reshape(encG, [R,C]);
enc_img(:,:,3) = reshape(encB, [R,C]);

% Stop timer
time_LR = toc;
fprintf('Lorenz–Rossler Encryption Time: %.6f seconds\n', time_LR);

%% Optional: Repeat multiple times for average time
num_runs = 100;
times_LR = zeros(1,num_runs);
for k = 1:num_runs
    tic;
    [X, Y, Z] = lorenz_rossler_system(delta, r, a, beta, b, c, [x0; y0; z0], N);
  Xn = uint8(mod(abs(X)*1e14, 256));
    Yn = uint8(mod(abs(Y)*1e14, 256));
Zn = uint8(mod(abs(Z)*1e14, 256));

% Extract and flatten image channels
Rvec = img(:,:,1);
Rvec = Rvec(:);

Gvec = img(:,:,2);
Gvec = Gvec(:);

Bvec = img(:,:,3);
Bvec = Bvec(:);

    [~, idxX] = sort(Xn); [~, idxY] = sort(Yn); [~, idxZ] = sort(Zn);
    shfR = Rvec(idxX); shfG = Gvec(idxY); shfB = Bvec(idxZ);
    encR = bitxor(shfR, Xn); encG = bitxor(shfG, Yn); encB = bitxor(shfB, Zn);
    enc_img = zeros(R,C,3,'uint8'); enc_img(:,:,1) = reshape(encR,[R,C]); enc_img(:,:,2) = reshape(encG,[R,C]); enc_img(:,:,3) = reshape(encB,[R,C]);
    times_LR(k) = toc;
end
avg_time_LR = mean(times_LR);
fprintf('Average Lorenz–Rossler Encryption Time over %d runs: %.6f seconds\n', num_runs, avg_time_LR);


