function decrypted_img = reverse_diffusion(encrypted_img, chaotic_seq)
%REVERSE_DIFFUSION Reverses the XOR-based diffusion of an image
%
% Usage:
%   decrypted_img = reverse_diffusion(encrypted_img, chaotic_seq)
%
% Inputs:
%   encrypted_img - Encrypted image after diffusion (RGB)
%   chaotic_seq   - Chaotic sequence used for diffusion (vector)
%
% Output:
%   decrypted_img - Image after reversing diffusion
%
% Description:
%   This function undoes the XOR-based diffusion operation performed
%   during encryption using the same chaotic sequence. The first pixel
%   uses an initialization vector (IV), and subsequent pixels are
%   recovered iteratively using the XOR relationship.

    [R,C,CH] = size(encrypted_img);
    N = R * C * CH;

    % Flatten image and chaotic key
    enc_vec = encrypted_img(:);
    K2 = chaotic_seq(:);

    % Compute initialization vector
    IV = uint8(mod(sum(K2), 256));

    % Preallocate reverse diffusion vector
    rev_diff = zeros(N,1,'uint8');

    % Reverse diffusion: first pixel
    rev_diff(1) = bitxor(enc_vec(1), bitxor(K2(1), IV));

    % Remaining pixels
    for i = 2:N
        rev_diff(i) = bitxor(enc_vec(i), bitxor(K2(i), enc_vec(i-1)));
    end

    % Reshape to image
    decrypted_img = reshape(rev_diff, [R, C, CH]);
end
