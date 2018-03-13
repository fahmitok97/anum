function y = forward_elim(L, b)
    [m,n] = size(L);
    % bandwidth dibawah diagonal utama
    q = bandwidth(L, 'lower');
    for j = 1 : n
        for i = 1 + j : min(n, j + q)
            b(i) = b(i) - L(i, j) * b(j);
        end
    end
    y = b;
end