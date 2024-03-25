clear;clc;
n = [4, 2]; L = 1; M = 12; k0a1 = [0.5 , 2.5]; k0di = [1.5]; thi = [0];

T = double(TmatrixL1TM(M, L, k0a1, k0di, thi, n));

fin = linspace(0, 2*pi, 60);

fi = linspace(0, 2*pi, 500);
sigb = zeros(size(fin));

for i = 1:length(fin)
      fn = fin(i);
      f0 = fn + pi;
      
    fop = 0;
    fop1 = 0;
    fop2 = 0;

    for n = -M:M
       for np = -M:M
          fop2 = fop2 + (-1j)^n * (1j)^np * exp(1j*(n * fn - np * fn)) * T(n + M + 1, np + M + 1);
          fop1 = fop1 + (-1j)^n * (1j)^np * exp(1j*(n * fi - np * fn)) * T(n + M + 1, np + M + 1); 
          
          fop = fop + (-1j)^n * (1j)^np * exp(1j*(n * f0 - np * fn)) * T(n + M + 1, np + M + 1); 
       end
    end
    
    fop = (1-1j)/sqrt(pi) * fop;
    sigb(i) = 2*pi*(abs(fop)^2);
    
    fop1 = (1-1j)/sqrt(pi) * fop1;
    fop2 = (1-1j)/sqrt(pi) * fop2;

    sigs = trapz(fi, ((abs(fop1)).^2));
    sige = -2 * sqrt(pi) * real((1+1j)*fop2);

    abs(sigs-sige)
end

plot(fin, sigb);