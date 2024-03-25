function dist = distij(I, J, k0di, thi)
    if(I == 0)
       dist = k0di(J);
       return;
    end
    if(J == 0)
       dist = k0di(I);
       return;
    end

    xi = k0di(I) * cos(thi(I));
    yi = k0di(I) * sin(thi(I));
    
    xj = k0di(J) * cos(thi(J));
    yj = k0di(J) * sin(thi(J));
    
    dist = sqrt((xi - xj)^2 + (yi - yj)^2);    
end