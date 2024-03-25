function ang = angij(I, J, k0di, thi)
    if (I == 0)
       ang = thi(J);
       return;
    end
    if (J == 0)
        ang = thi(I) + pi;
        if thi(I) > pi
            ang = thi(I) - pi;
        end
        return;
    end
    
    
    xi = k0di(I) * cos(thi(I));
    yi = k0di(I) * sin(thi(I));
    
    xj = k0di(J) * cos(thi(J));
    yj = k0di(J) * sin(thi(J));
    
    ang = atan2(yi - yj, xi - xj);
    
    if ang < 0
        ang = ang + 2*pi;
    end
    
end