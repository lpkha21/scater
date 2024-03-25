function H = besH(M, I, J, k0ai, ni)
    if I == 0
       rad = k0ai(J); 
    else
        rad = ni(I) * k0ai(J);
    end
    
    H = zeros(2*M + 1);
    
    for m = -M:M
       H(m + M + 1, m + M + 1) = besselh(m, 1, rad); 
    end

end