function decrypted_img = reverse_confusion(confused_img, chaotic_seq)
%REVERSE_CONFUSION Reverses the pixel permutation (confusion) of an image
%
% Usage:
%   decrypted_img = reverse_confusion(confused_img, chaotic_seq)
%
% Inputs:
%   confused_img  - Image after confusion (RGB)
%   chaotic_seq   - Chaotic sequence used for permutation (vector)
%
% Output:
%   decrypted_img - Image after reversing confusion
%
% Description:
%   This function restores the original pixel positions of an image that
%   has been scrambled using a chaotic sequence. The inverse permutation
%   is computed using the same chaotic sequence that was used during
%   encryption.

    [R,C,CH] = size(confused_img);
    N = R * C * CH;

    % Flatten chaotic sequence and compute permutation indices
    K1 = double(chaotic_seq(:));
    [~, perm_index] = sort(K1);

    % Apply inverse permutation
    dec_vec = zeros(N,1,'uint8');
    dec_vec(perm_index) = confused_img(:);

    % Reshape to image
    decrypted_img = reshape(dec_vec, [R, C, CH]);
end
