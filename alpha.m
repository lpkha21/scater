function a = alpha(M, I, J, k0di, thi)

    [n, m] = meshgrid(-M:M, -M:M);
    
    if I < J
        dist = distij(I, J, k0di, thi);
        ang = angij(I, J, k0di, thi);

        a = besselh(n-m, 1, dist) .* exp(-1j * (n - m) * ang);
    
    else
        dist = distij(J, I, k0di, thi);
        ang = angij(J, I, k0di, thi);
        
        a = (-1) .^ (n-m) .* besselh(n - m, 1, dist) .* exp(-1j * (n-m) * ang);
    end
end