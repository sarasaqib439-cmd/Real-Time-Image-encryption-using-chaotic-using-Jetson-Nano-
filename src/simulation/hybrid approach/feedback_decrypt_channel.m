function P = feedback_decrypt_channel(Cc, Kk)
    n = numel(Cc);
    P = zeros(n,1,'uint8');
    IV = uint8(Kk(1));
    P(1) = bitxor(uint8(Cc(1)), bitxor(uint8(Kk(1)), uint8(IV)));
    for jj = 2:n
        % P(i) = C(i) XOR K(i) XOR C(i-1)
        P(jj) = bitxor(uint8(Cc(jj)), bitxor(uint8(Kk(jj)), uint8(Cc(jj-1))));
    end
end
