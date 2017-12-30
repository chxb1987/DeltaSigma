N = 32768; OSR = N/66; M = 3; nlev = M+1;
H = synthesizeNTF(2,OSR,1);
fB = ceil(N/(2*OSR)); ftest = floor(2/3*fB);
u = 1.5*sin(2*pi*ftest/N*(0:N-1));	% half-scale sine-wave input
[v_tmp,xn,xnmax,y] = simulateDSM(u,H,nlev); 
f = linspace(0,0.5,N/2+1);
sigma_d = 0.01;		% 1% mismatch
window = ds_hann(N)/(M*N/8);
windowM = repmat(window,M,1);
spec = fft(v_tmp.*window);

mtf = synthesizeNTF(2,OSR,1,2);	    % Second-order shaping
fin = round(0.01*N);
if isempty(mtf)
    sv = ds_therm(v_tmp,M);
else
    sv = simulateMS(v_tmp,M,mtf,0);
end
Sdd = sum(abs(fft(sv.*windowM,[],2)).^2);

figure;
plot(log10(f),dbv(spec(1:N/2+1))); hold on;
snr = calculateSNR(spec(3:fB+1),ftest-2);
plot(log10(f),dbv(sqrt(Sdd(1:N/2+1))));
snr_d = calculateSNR(sqrt(Sdd(3:fB+1)),ftest-2);
NBW = 1.5/N;
Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / 3;
plot(log10(f),dbp(Sqq*NBW),'m','Linewidth',2); grid on; 
text(0.001, -90, sprintf('NBW = %4.1E x f_s ',NBW),'Hor','right');
legend(sprintf('2-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),sprintf('2-bit Noise Shaped SNR = %4.1fdB @ OSR = %d',snr_d,OSR),'2-bit Expected');

          



