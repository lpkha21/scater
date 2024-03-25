clear;clc;
ep0 = 1/36/pi * 1e-9; mi0 = 4*pi*1e-7;
ni = [4,4];
ai = [0.5 2.5];
% ni0 = sqrt(mi0/ep0); ni1 = ni0/4; ni2 = ni0/2;

ni0 = sqrt(ep0/mi0); ni1 = ni0*4; ni2 = ni0*2;
first = ni2*dbesH(12,2,1,ai,ni) - ni1*dbesJ(12,1,1,ai,ni)/besJ(12,1,1,ai,ni)*besH(12,2,1,ai,ni);
second = ni2*dbesJ(12,2,1,ai,ni) - ni1*dbesJ(12,1,1,ai,ni)/besJ(12,1,1,ai,ni)*besJ(12,2,1,ai,ni);

T11 = -first\second;

[n, m] = meshgrid(-12:12, -12:12);

b01 = besselj(n - m, 1.5);
b10 = besselj(n-m, 1.5) * exp(-1j *(n-m)*pi);

t1 = b01 * T11 / b01;

f = besJ(12,2,2,ai,ni)  + besH(12,2,2,ai,ni) * t1;
df = dbesJ(12,2,2,ai,ni)  + dbesH(12,2,2,ai,ni) * t1;

TLfirst = ni0*dbesH(12,0,2,ai,ni) - ni2*df/f * besH(12,0,2,ai,ni);
TLsecond = ni0*dbesJ(12,0,2,ai,ni) - ni2*df/f * besJ(12,0,2,ai,ni);

TL1 = -(TLfirst\TLsecond);

s = 0;
fi = linspace(0, 2*pi, 500);

fs = 0;
fe = 0;
s0 = 0;
fin = linspace(0, 2*pi);
for n = -12:12
   for np = -12:12
        fs = fs + (-1j)^n * (1j)^np * exp(1j*(n*fi - np*pi)) * TL1(n + 13, np + 13); 
        fe = fe + (-1j)^n * (1j)^np * exp(1j*(n*pi - np*pi)) * TL1(n + 13, np + 13); 
        s = s + (-1j)^n * (1j)^np * exp(1j*pi*(2*n - np)) * TL1(n + 13, np + 13); 
        s0 = s0 + (-1j)^n * (1j)^np * exp(1j*((fin+pi)*n - fin*np)) * TL1(n + 13, np + 13); 
   end
end
fs = (1-1j)/sqrt(pi) * fs;
fe = (1-1j)/sqrt(pi) * fe;

sigs = trapz(fi, abs(fs).^2);
sige = -2*sqrt(pi)*real((1+1j)*fe);
disp("Optical Theorem")
abs(sigs-sige)

disp("k0*sigb at fin = pi")
2 * abs((1-1j)*s)^2

plot(fin, 2 * abs((1-1j)*s0).^2);xlim([0 2*pi]);


