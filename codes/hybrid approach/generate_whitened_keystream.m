function ks = generate_whitened_keystream(source_vec, N_needed, blk_d)
    % source_vec: normalized double vector (length L)
    % produces uint8 keystream of length N_needed using SHA-256 per block
    ks = zeros(N_needed,1,'uint8');
    pos = 1;
    idx = 1;
    L = length(source_vec);
    while pos <= N_needed
        i1 = ( (idx-1)*blk_d + 1 );
        blockIdx = mod((i1-1):(i1+blk_d-2), L) + 1;
        block = double(source_vec(blockIdx));
        hbytes = sha256_bytes_of_doubles(block); % 32 bytes
        take = min(length(hbytes), N_needed - pos + 1);
        ks(pos:pos+take-1) = hbytes(1:take);
        pos = pos + take;
        idx = idx + 1;
    end
end

