function Tl = TiLbTM(M, L, l,ind, k0di, k0ai, thi, ni)
    if l == 1
        Tl = TiLbTMbasecase(M,L,  ind, k0ai, ni) * beta(M, ind, 0, k0di, thi);
        return;
    end

    first_sum = 0;
    
    for i = 1:l-1
        first_sum = first_sum + alpha(M, l, i, k0di, thi) * TiLbTM(M, L, l - 1, i, k0di, k0ai,  thi, ni) * ...
            beta(M, 0, i, k0di, thi) * alpha(M, i, l, k0di, thi);
    end

    first_part = eye(2*M + 1) - TiLbTMbasecase(M,L,  l, k0ai, ni) * first_sum;
    first_part = first_part \  TiLbTMbasecase(M,L,  l, k0ai, ni);
    
    second_sum = 0;
    
    for i = 1:l-1
        second_sum = second_sum + alpha(M, l, i, k0di, thi) * TiLbTM(M, L, l-1, i, k0di, k0ai, thi, ni);
    end
    
    second_part = beta(M, l, 0, k0di, thi) + second_sum;
    
    TLB = first_part * second_part;
    
    if (ind == l)
        Tl = TLB;
       return; 
    end
    
    
    Tl = TiLbTM(M, L, l-1, ind, k0di, k0ai, thi, ni)* (eye(2*M + 1) + beta(M, 0, ind, k0di, thi) * alpha(M, ind, l, k0di, thi) * TLB);
end