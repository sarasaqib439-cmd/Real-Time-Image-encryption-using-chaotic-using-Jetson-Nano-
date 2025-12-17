function diffused_img = diffusion_stage(confused_img, chaotic_seq)
%DIFFUSION_STAGE Performs XOR-based diffusion on a confused image
%
% Inputs:
%   confused_img - Image after confusion
%   chaotic_seq  - Chaotic sequence (length = numel(confused_img))
%
% Output:
%   diffused_img - Encrypted image after diffusion

    [R, C, CH] = size(confused_img);
    N = R * C * CH;

    confused_vec = confused_img(:);
    K2 = chaotic_seq(:);

    diff_vec = zeros(N,1,'uint8');

    % Initialization vector
    IV = uint8(mod(sum(K2), 256));

    % First pixel
    diff_vec(1) = bitxor(confused_vec(1), bitxor(K2(1), IV));

    % Remaining pixels
    for i = 2:N
        diff_vec(i) = bitxor(confused_vec(i), bitxor(K2(i), diff_vec(i-1)));
    end

    % Reshape back to image
    diffused_img = reshape(diff_vec, [R C CH]);
end
