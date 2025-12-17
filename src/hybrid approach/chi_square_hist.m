function chi2vals = chi_square_hist(img, method_name)
    % img: input RGB image (uint8)
    % method_name: string specifying the encryption method (e.g., 'Lorenz–Rossler')
    % chi2vals: 1x3 array containing chi-square values for R, G, B channels
    
    if nargin < 2
        method_name = 'Image'; % default name if not provided
    end
    
    chi2vals = zeros(1,3);  % preallocate
    
    for ch = 1:3
        channel_data = img(:,:,ch);
        channel_data = channel_data(:);       % flatten
        N = numel(channel_data);              % total number of pixels
        counts = histcounts(channel_data, 0:256); % histogram counts
        E = N / 256;                          % expected count per bin
        chi2vals(ch) = sum((counts - E).^2 / E); % chi-square formula
    end
    
    % Display results
    fprintf('Chi-Square Values (%s):\n', method_name);
    fprintf('Red Channel   : %.4f\n', chi2vals(1));
    fprintf('Green Channel : %.4f\n', chi2vals(2));
    fprintf('Blue Channel  : %.4f\n', chi2vals(3));
end
