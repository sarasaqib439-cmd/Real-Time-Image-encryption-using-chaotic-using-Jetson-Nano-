function confused_img = confusion_stage(img, chaotic_seq)
%CONFUSION_STAGE Scrambles an image based on a chaotic sequence
%
% Inputs:
%   img         - Input image (RGB)
%   chaotic_seq - Chaotic sequence (length = numel(img))
%
% Output:
%   confused_img - Pixel-permuted (confused) image

    [R, C, CH] = size(img);
    N = R * C * CH;

    img_vec = img(:);
    K1 = double(chaotic_seq(:));

    % Sort chaotic sequence to get permutation indices
    [~, perm_index] = sort(K1);

    % Apply permutation
    confused_vec = img_vec(perm_index);

    % Reshape back to image
    confused_img = reshape(confused_vec, [R C CH]);
end
