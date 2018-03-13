function x = backwardSub(U, b)
    [m,n] = size(U);
    % bandwidth diatas diagonal utama
    p = bandwidth(U, 'upper');
    for j = n : -1 : 1
        b(j) = b(j) / U(j, j);
        for i = max(1, j - p) : j - 1
            b(i) = b(i) - U(i, j) * b(j);
        end
    end
    x = b;
end