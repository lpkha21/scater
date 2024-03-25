clear;clc;
ka = 2.5; n1 = 2; fin = 3*pi/2;
ep0 = 1/36/pi * 1e-9; mi0 = 4*pi*1e-7; ni0 = sqrt(ep0/mi0); ni1 = n1*ni0;


fi0 = linspace(0, 2*pi);

f = besJ(12,1,1,ka,n1) ;
df = dbesJ(12,1,1,ka,n1) ;

TLfirst = ni0*dbesH(12,0,1,ka,n1) - ni1*df/f * besH(12,0,1,ka,n1);
TLsecond = ni0*dbesJ(12,0,1,ka,n1) - ni1*df/f * besJ(12,0,1,ka,n1);

TL1 = -(TLfirst\TLsecond);

s = 0;
for n = -12:12
   for np = -12:12
        s = s + (-1j)^n * (1j)^np * exp(1j*(fi0*n - fin*np)) * TL1(n + 13, np + 13); 
   end
end

sig = 2*abs((1-1j)*s).^2;

plot(fi0, sig)
xlim([0 2*pi]); ylim([0 140])
