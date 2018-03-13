function [L,U,P] = luFactPivot(A)
% LU factorization with partial pivoting

    [n,n] = size(A);
    [p, q] = bandwidth(A);
    L = eye(n);
    P = L;
    U = A;
    for k = 1:n-1
        [pivot m] = max(abs(U(k:min(k+p, n), k))); % elemen maks kolom k
        m = m+k-1;
        if m~=k
            % tukar row m and k di U
            temp = U(k,:);
            U(k,:) = U(m,:);
            U(m,:) = temp;

            % tukar row m and k di P
            temp = P(k,:);
            P(k,:) = P(m,:);
            P(m,:) = temp;
            if k >= 2
                temp = L(k,1:k-1);
                L(k,1:k-1) = L(m,1:k-1);
                L(m,1:k-1) = temp;
            end
        end

        % calculate LU bandwidth
        for j=k+1:min(k+p, n)
            L(j,k) = U(j,k)/U(k,k);
            U(j,k+1:n) = U(j,k+1:n)-L(j,k)*U(k,k+1:n);
        end
    end
end