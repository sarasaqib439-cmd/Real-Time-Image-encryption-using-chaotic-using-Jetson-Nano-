function C = feedback_encrypt_channel(P, K)
    % P: plaintext column vector (uint8 or numeric)
    % K: keystream uint8 vector (same length)
    n = numel(P);
    C = zeros(n,1,'uint8');
    IV = uint8(K(1)); % use first keystream byte as IV
    C(1) = bitxor(uint8(P(1)), bitxor(uint8(K(1)), uint8(IV)));
    for ii = 2:n
        % CBC-like: C(i) = P(i) XOR K(i) XOR C(i-1)
        C(ii) = bitxor(uint8(P(ii)), bitxor(uint8(K(ii)), C(ii-1)));
    end
end

