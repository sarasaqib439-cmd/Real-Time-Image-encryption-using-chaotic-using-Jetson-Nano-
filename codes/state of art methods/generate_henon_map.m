function henon_seq = generate_henon_map(a, b, x0, y0, N)
%GENERATE_HENON_MAP Generates a Henon chaotic sequence
%
% Inputs:
%   a, b  - Henon map parameters (e.g., a=1.4, b=0.3)
%   x0    - initial x value (0 < x0 < 1)
%   y0    - initial y value (0 < y0 < 1)
%   N     - length of the sequence (number of pixels)
%
% Output:
%   henon_seq - normalized chaotic sequence (uint8, range 0-255)

    x = zeros(1, N);
    y = zeros(1, N);

    x(1) = x0;
    y(1) = y0;

    for i = 2:N
        x(i) = 1 - a * x(i-1)^2 + y(i-1);
        y(i) = b * x(i-1);
    end

    % Normalize to [0, 255] as uint8
    henon_seq = uint8(255 * mat2gray(x));

end
