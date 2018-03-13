function [x] = solveLUFact(A, b)
    %% Expecting 2 arguments, a square matrix and a vector
    if (nargin ~= 2)
        print_usage ();
    end

    % check if A is a nxn matrix
    if (!issquare(A))
        error('Only square systems')
    end

    % check if A is a banded matrix
    % its upper and lower must not equal to its size - 1
    [n,n] = size(A);
    [upper, lower] = bandwidth(A);
    if (upper == n-1 && lower == n-1)
        error('A has to be a banded nxn matrix')
    end
    % check if B is in the same dimension as A
    if (size(A,1) != size(b,1))
            error('b must be in the same dimension as A')
    end

    [L, U, P] = luFactPivot(A) % LU factorization of A
    y = forwardElim(L, P*b) % solve Ly = Pb;
    x = backwardSub(U, y) % solve Ux = y;
end