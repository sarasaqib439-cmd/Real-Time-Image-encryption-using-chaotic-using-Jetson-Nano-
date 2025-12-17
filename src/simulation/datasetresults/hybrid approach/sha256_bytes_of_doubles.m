function outBytes = sha256_bytes_of_doubles(dvals)
    % dvals: double vector -> return uint8(32) digest
    db = typecast(dvals(:),'uint8');         % 8 bytes per double
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(db);
    digest = md.digest();                    % Java byte[]
    % Convert Java byte[] to MATLAB uint8 (handles sign)
    outBytes = uint8(zeros(1, numel(digest)));
    for k = 1:numel(digest)
        % Java bytes are signed -128..127, convert to 0..255
        outBytes(k) = uint8(int16(digest(k)) + 256 * (digest(k) < 0));
    end
end