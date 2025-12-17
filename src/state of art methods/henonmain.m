%% Chaotic Maps Image Encryption
clc; clear; close all;

%% HENON CHAOTIC MAP GENERATION
img =imread(' C:\Users\rimsha malik\Desktop\PhD NUST\computer vision\assignment 2\project\codes\state of art methods\paper.PNG');
a = 1.4;   
b = 0.3;  
x0 = 0.9; 
y0 = 0.5; 

[R, C, CH] = size(img);
N = R*C*CH;


%% Encryption Process:
% Generate Henon chaotic sequence
henon_seq = generate_henon_map(a, b, x0, y0, N);


%% CONFUSION USING HENON
% Confusion
confused_henon = confusion_stage(img, henon_seq );

%% DIFFUSION USING HENON
% Diffusion
encrypted_henon = diffusion_stage(confused_henon , henon_seq);
figure;
subplot(131), imshow(img), title('Original');
subplot(132), imshow(confused_henon), title('Confused Image (Henon)');
subplot(133), imshow(encrypted_henon), title('Encrypted Image (Henon)');

%% %% Histogram Analysis for Henon Map Encryption
fprintf('=== Henon Map Encryption Histogram Analysis ===\n');
plot_hist(img, 'Original Image');
plot_hist(confused_henon, 'Confused Image (Henon)');
plot_hist(encrypted_henon, 'Encrypted Image (Henon)');
chi_Image   = chi_square_hist(img, 'Original Image');
chi_henon    = chi_square_hist(encrypted_henon,  'Encrypted Image (Henon)');
%% Decryption process:
a = 1.4;   
b = 0.3;  
x0 = 0.9; 
y0 = 0.5; 


% Generate Henon chaotic sequence
henon_seq = generate_henon_map(a, b, x0, y0, N);
%% 
% Step 1: Reverse diffusion
after_diffusion = reverse_diffusion(encrypted_henon, henon_seq );

% Step 2: Reverse confusion
decrypted_henon = reverse_confusion(after_diffusion, henon_seq );

% Display
figure;
subplot(131), imshow(encrypted_henon), title('Encrypted Image');
subplot(132), imshow(after_diffusion), title('After Reverse Diffusion');
subplot(133), imshow(decrypted_henon), title('Decrypted Image');

%% Correlation Analysis:
correlation_analysis(img, 'Original Image Correlation');
correlation_analysis(encrypted_henon, 'Henon Encrypted Image Correlation');

%% %% Sensitivity of Key for Decryption process:

%% For Henon:
% Encryption Key
a = 1.4;   
b = 0.3;  
x0 = 0.9; 
y0 = 0.5; 

[R, C, CH] = size(img);
N = R*C*CH;

% Generate Henon chaotic sequence
henon_seq = generate_henon_map(a, b, x0, y0, N);
%% CONFUSION USING HENON
confused_henon = confusion_stage(img, henon_seq );

%% DIFFUSION USING HENON
encrypted_henon = diffusion_stage(confused_henon , henon_seq);
figure;
subplot(131), imshow(img), title('Original');
subplot(132), imshow(encrypted_henon), title('Encrypted Image (Henon)');

%%  Decryption Key
a = 1.4;   
b = 0.3;  
x0 = 0.9; 
y0 = 0.50000001; 


% Generate Henon chaotic sequence
henon_seq = generate_henon_map(a, b, x0, y0, N);


% Step 1: Reverse diffusion
after_diffusion = reverse_diffusion(encrypted_henon, henon_seq );

% Step 2: Reverse confusion
decrypted_henon = reverse_confusion(after_diffusion, henon_seq );

% Display
subplot(133), imshow(decrypted_henon), title('Decrypted Image');
sgtitle(' Image encrypted key:(0.9, 0.5) & Decrypted image key: (0.9, 0.50000001)')

%% %% Timing Analysis for Henon Map Encryption
fprintf('=== Henon Map Encryption Timing Analysis ===\n');

% Start timer
tic;

% Generate Henon chaotic sequence
henon_seq = generate_henon_map(a, b, x0, y0, N);

% Confusion
confused_henon = confusion_stage(img, henon_seq );

% Diffusion
encrypted_henon = diffusion_stage(confused_henon , henon_seq);

% Stop timer
time_henon = toc;
fprintf('Henon Map Encryption Time: %.6f seconds\n', time_henon);

%% Optional: Repeat multiple times for average time
num_runs = 100;
times_henon = zeros(1,num_runs);
for k = 1:num_runs
    tic;
    henon_seq = generate_henon_map(a, b, x0, y0, N);
    confused_henon = confusion_stage(img, henon_seq );
    encrypted_henon = diffusion_stage(confused_henon , henon_seq);
    times_henon(k) = toc;
end
avg_time_henon = mean(times_henon);
fprintf('Average Henon Map Encryption Time over %d runs: %.6f seconds\n', num_runs, avg_time_henon);

