function dJ = dbesJ(M, I, J, k0ai, ni)
    if I == 0
       rad = k0ai(J); 
    else
        rad = ni(I) * k0ai(J);
    end
    
    dJ = zeros(2*M + 1);
    
    for m = -M:M
       dJ(m + M + 1, m + M + 1) = m*besselj(m, rad)/rad - besselj(m + 1, rad); 
    end

end