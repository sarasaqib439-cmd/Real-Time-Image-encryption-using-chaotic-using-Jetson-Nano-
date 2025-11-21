function tent_seq = generate_tent_map(alpha, x0, N)
%GENERATE_TENT_MAP Generates a Tent chaotic sequence
%
% Inputs:
%   alpha - control parameter of Tent map (e.g., 1.99)
%   x0    - initial condition / secret key (0 < x0 < 1)
%   N     - length of the sequence (number of pixels)
%
% Output:
%   tent_seq - normalized chaotic sequence (uint8, range 0-255)

    tent = zeros(1, N);
    tent(1) = x0;

    for i = 2:N
        if tent(i-1) < 0.5
            tent(i) = alpha * tent(i-1);
        else
            tent(i) = alpha * (1 - tent(i-1));
        end
    end

    % Normalize to [0, 255] as uint8
    tent_seq = uint8(255 * mat2gray(tent));

end
