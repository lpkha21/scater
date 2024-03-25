function Tl = TmatrixL1TM(M, L, k0ai, k0di, thi, ni)
ep0 = 8.8541878128 * 1e-12;
mi0 = 1.25663706212 * 1e-6;
niu0 = vpa(sqrt(ep0/mi0),100);
niuL = vpa(ni(L+1)*niu0,100);

Tm = 0;

for i = 1:L
   Tm = Tm + vpa(beta(M, 0, i, k0di, thi), 100) * vpa(TiLbTM(M, L, L, i, k0di, k0ai, thi, ni), 100);
end

f = vpa(besJ(M,L+1,L+1,k0ai,ni),100) + vpa(besH(M,L+1,L+1,k0ai,ni) * Tm,100);
df = vpa(dbesJ(M,L+1,L+1,k0ai,ni),100) + vpa(dbesH(M,L+1,L+1,k0ai,ni) * Tm,100);

first  = niu0*dbesH(M,0,L+1,k0ai,ni) - niuL*((vpa(df,100) / vpa(f,100)) * besH(M,0,L+1,k0ai,ni));
second = niu0*dbesJ(M,0,L+1,k0ai,ni) - niuL*((vpa(df,100) / vpa(f,100)) * besJ(M,0,L+1,k0ai,ni));

Tl = -(vpa(first,100) \ vpa(second,100));
end