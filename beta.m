function b = beta(M, I, J, k0di, thi)    
    
    [n, m] = meshgrid(-M:M, -M:M);
    
    if I < J
      dist = distij(I, J, k0di, thi);
      ang = angij(I, J, k0di, thi);

      b = besselj(n-m, dist) .* exp(-1j * (n - m) * ang);
    
    else
        dist = distij(J, I, k0di, thi);
        ang = angij(J, I, k0di, thi);
        
        b = (-1) .^ (n-m) .* besselj(n - m, dist) .* exp(-1j * (n-m) * ang);
    end
    
end