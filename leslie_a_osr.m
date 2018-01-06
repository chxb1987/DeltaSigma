m = 31:2:37; snr_mx = zeros(1,numel(m)); osr_mx = zeros(1,numel(m));
for k=1:numel(m)
N = 2.^(7:20); nlev = 2; snrmax = zeros(1,numel(N));
for j=1:numel(N)
OSR = N(j)/(2*m(k));
H = synthesizeNTF(2,OSR,1);
fB = ceil(N(j)/(2*OSR)); ftest = floor(2/3*fB);
u = 0.5*sin(2*pi*ftest/N(j)*(0:N(j)-1));	% half-scale sine-wave input
[v_tmp,xn,xnmax,y] = simulateDSM(u,H,nlev); 
f = linspace(0,0.5,N(j)/2+1);

v_tmp_a = v_tmp./2;

snr_lsli = zeros(1,19);
for i=2:20
v = ds_quantize(((2^i)-1).*y,2^i); 
v_a = v./(2^i);

v_tmp_1_a = 0.775.*[0 v_tmp_a(1:numel(v)-1)] - 0.5585.*[0 0 v_tmp_a(1:numel(v)-2)];

v_1_a = 0.5.*(v_a - [0 v_a(1:numel(v)-1)]);
v_2_a = 0.5.*(v_1_a - [0 v_1_a(1:numel(v)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;

spec = fft(v_tmp.*ds_hann(N(j)))/(N(j)/4);
% figure;
% subplot(2,1,1);
% plot(log10(f),dbv(spec(1:N(j)/2+1))); hold on;
snr = calculateSNR(spec(3:fB+1),ftest-2);
NBW = 1.5/N(j);
Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / 3;
% plot(log10(f),dbp(Sqq*NBW),'m','Linewidth',2); grid on;
% text(0.001, -90, sprintf('NBW = %4.1E x f_s ',NBW),'Hor','right');
% legend(sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),'1-bit Expected');

% subplot(2,1,2);
spec_v = fft(v_a.*ds_hann(N(j)))/(N(j)/4);
% plot(log10(f),dbv(spec_v(1:N(j)/2+1))); hold on;
snr_v = calculateSNR(spec_v(3:fB+1),ftest-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N(j)))/(N(j)/4);
% plot(log10(f),dbv(spec(1:N(j)/2+1))); hold on;
% plot(log10(f),dbv(spec_lsli(1:N(j)/2+1))); grid on;
snr_lsli(i-1) = calculateSNR(spec_lsli(3:fB+1),ftest-2);

% legend(sprintf('2-bit Quantizer SNR = %4.1fdB @ OSR = %d',snr_v,OSR),sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),sprintf('Leslie Singh SNR = %4.1fdB @ OSR = %d',snr_lsli,OSR));
end
snrmax(j) = max(snr_lsli-snr);  
end
% plot(N./22,snrmax); grid on;
snr_mx(k) = max(snrmax); osr_mx(k) = N(find(snrmax == max(snrmax)))/(2*m(k));
end
figure;
plot(osr_mx,snr_mx); grid on;