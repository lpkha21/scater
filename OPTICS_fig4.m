clear;clc;
L = 0; ka = 2*pi*0.03; n = 3;

M = 12;
t = double(TmatrixL1TM(M ,0,ka, 0, 0, n));

fi = pi;
f0 = linspace(0, 2*pi, 40);
sigd = zeros(size(f0));
for v = 1:length(f0)
    fop = 0;
    for n = -M:M
       for np = -M:M
          fop = fop + (-1j)^n * (1j)^np * exp(1j*(n * f0(v) - np * fi)) * t(n + M + 1, np + M + 1); 
       end
    end

   fop = (1-1j)/sqrt(pi) * fop;

   sigd(v) = 2*pi*abs(fop)^2;
end


plot(f0, db(sigd))

ylim([-30, 40])