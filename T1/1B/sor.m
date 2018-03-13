function [x,iter]=sor(A,b)
    %% Expecting 2 arguments, a square matrix and a vector
    if nargin ~= 2
        print_usage();
    end

    [n,m] = size(A);
    if n ~= m
        error('Only square systems');
    end

    [k] = size(b);
    if k ~= m
        error('Size of b is not equal to number of column of A');
    end

    %% CONSTANT BEGIN %%
    MAX_ITERATION = 10000;
    EPS = 1e-9;
    INITIAL_ERR=100;
    %% CONSTANT END %%

    iter = 0;
    x = ones(n,1);
    omega = omegaOptimal(m, n);
    err = INITIAL_ERR;

    while err > EPS && iter < MAX_ITERATION
        iter = iter + 1;

        for i = 1:n
            prev_x(i) = x(i);
        end

        for i = 1:n
            sigma = 0;
            for j = 1:n
                if i ~= j, sigma = sigma + A(i,j) * prev_x(j); end
            end
            x(i) = (omega * (b(i) - sigma)/A(i,i)) + ((1 - omega) * prev_x(i));
        end

        err = 0;
        for i = 1:n
            err = max(err, abs((x(i) - prev_x(i)) / x(i)) * 100);
        end
    end
    return
