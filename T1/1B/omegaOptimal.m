function w = omegaOptimal(m, n)
    c = cos(pi/m) + cos(pi/n);
    w = 4 / (2 + sqrt(4 + c*c));
end