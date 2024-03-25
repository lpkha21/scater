function J = besJ(M, I, J, k0ai, ni)
    if I == 0
       rad = k0ai(J); 
    else
        rad = k0ai(J) * ni(I);
    end
    
    J = zeros(2*M + 1);
    
    for m = -M:M
       J(m + M + 1, m + M + 1) = besselj(m, rad); 
    end

end