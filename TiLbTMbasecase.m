function base = TiLbTMbasecase(M,L,  ind, k0ai, ni)
    mi0 = 4*pi*1e-7;
    ep0 = 8.8541878128e-12;

    ni0 = sqrt(ep0/mi0);
    niL = ni(L+1) * ni0;
    niInd = ni(ind) * ni0;
    
    base1 = -(niL * dbesH(M, L + 1, ind, k0ai, ni) - niInd* (dbesJ(M, ind, ind, k0ai, ni)) / dbesJ(M, ind, ind, k0ai, ni) * besH(M, L + 1, ind, k0ai, ni) );
    
    base2 = niL * dbesJ(M, L + 1, ind, k0ai, ni) - niInd* (dbesJ(M, ind, ind, k0ai, ni)) / dbesJ(M, ind, ind, k0ai, ni) * besJ(M, L + 1, ind, k0ai, ni);
    
    base = base1 \ base2;    
end