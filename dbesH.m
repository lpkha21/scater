function dH = dbesH(M, I, J, k0ai, ni)
    if I == 0
       rad = k0ai(J); 
    else
        rad = ni(I) * k0ai(J);
    end
    
    dH = zeros(2*M + 1);
    
    for m = -M:M
       dH(m + M + 1, m + M + 1) = m*besselh(m, 1, rad)/rad - besselh(m + 1, 1, rad); 
    end

end