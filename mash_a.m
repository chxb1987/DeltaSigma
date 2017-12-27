N = 32768; OSR = N/66; nlev = 2;
H = synthesizeNTF(2,OSR,1);
fB = ceil(N/(2*OSR)); ftest = floor(2/3*fB);
u = 0.5*sin(2*pi*ftest/N*(0:N-1));	% half-scale sine-wave input
[v_tmp,xn,xnmax,y] = simulateDSM(u,H,nlev); 
f = linspace(0,0.5,N/2+1);

v_tmp_a = v_tmp./2;

snr_mash = zeros(1,19);
for i=2:20
[v,xn_v,xnmax_v,y_v] = simulateDSM(((2^i)-1).*(v_tmp-y),H,2^i);
v_a = v./(2^i);

v_tmp_1_a = 0.5.*(v_tmp_a - 1.225.*[0 v_tmp_a(1:numel(v)-1)] + 0.4415.*[0 0 v_tmp_a(1:numel(v)-2)]);

v_1_a = 0.25.*(v_a - [0 v_a(1:numel(v)-1)]);
v_2_a = 0.5.*(v_1_a - [0 v_1_a(1:numel(v)-1)]);

v_mash_a = v_tmp_1_a - v_2_a;

spec = fft(v_tmp.*ds_hann(N))/(N/4);
% figure;
% subplot(2,1,1);
% plot(log10(f),dbv(spec(1:N/2+1))); hold on;
snr = calculateSNR(spec(3:fB+1),ftest-2);
NBW = 1.5/N;
Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / 3;
% plot(log10(f),dbp(Sqq*NBW),'m','Linewidth',2); grid on;
% text(0.001, -90, sprintf('NBW = %4.1E x f_s ',NBW),'Hor','right');
% legend(sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),'1-bit Expected');

% subplot(2,1,2);
spec_v = fft(v_a.*ds_hann(N))/(N/4);
% plot(log10(f),dbv(spec_v(1:N/2+1))); hold on;
snr_v = calculateSNR(spec_v(3:fB+1),ftest-2);

spec_mash = fft(v_mash_a.*ds_hann(N))/(N/4);
% plot(log10(f),dbv(spec(1:N/2+1))); hold on;
% plot(log10(f),dbv(spec_mash(1:N/2+1))); grid on;
snr_mash(i-1) = calculateSNR(spec_mash(3:fB+1),ftest-2);

% legend(sprintf('2-bit Quantizer SNR = %4.1fdB @ OSR = %d',snr_v,OSR),sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),sprintf('MASH SNR = %4.1fdB @ OSR = %d',snr_mash(i-1),OSR));
end
plot(2:20,snr_mash-snr); grid on; 