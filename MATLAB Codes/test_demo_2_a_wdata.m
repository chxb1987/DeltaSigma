function [] = test_demo_2_a_wdata(optim)
N = 32768; OSR = N/66; nlev = 2;
fB = ceil(N/(2*OSR)); ftest = floor(2/3*fB);
u = 0.7071*sin(2*pi*fB/N*(0:N-1));	% half-scale sine-wave input
sdin = csvread('Q_w2bit.csv',1,0);
sdin(:,2) = 2.*(double(sdin(:,2) > 0.9))-1;
v_tmp = sdin(10:(N+9),2)';

% optim = [1,1.5,0.75];
z = zpk('z',1);
H = ((z-optim(1))^2)/(z^2-optim(2)*z+optim(3));
% [~,~,~,y] = simulateDSM(u,H,nlev);

yin1 = csvread('von_w2bit.csv',1,0);
% yin2 = csvread('vip1_op3.csv',1,0);
yin1(:,2) = (yin1(:,2) - 0.9)./0.1;
% yin2(:,2) = (yin2(:,2) - 0.9)./0.125;
y = yin1(10:(N+9),2)';
% [v_tmp,xn,xnmax,y] = simulateDSM(u,H,nlev); 
f = linspace(0,0.5,N/2+1);

% v_tmp_a = (v_tmp - 0.9)./0.9;
% y = 5.*(y - 0.9);
v_tmp_a = v_tmp./2;

snr_lsli = zeros(1,19);
for i=2:20
v = ds_quantize(((2^i)-1).*y,2^i); 
v_a = v./((2^i));
% v_a = v_a./16;

v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 v_tmp_a(1:numel(v)-1)] - (optim(1)^2-optim(3)).*[0 0 v_tmp_a(1:numel(v)-2)];

v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v)-1)]);

% v_1_a = v_a - 2.0554.*[0 v_a(1:numel(v)-1)] + 1.0562.*[0 0 v_a(1:numel(v)-2)];

v_lsli_a = v_tmp_1_a + v_2_a;

spec = fft(v_tmp_a.*ds_hann(N))/(N/4);
% figure;
% subplot(2,1,1);
% plot(log10(f),dbv(spec(1:N/2+1))); hold on; grid on;
snr = calculateSNR(spec(3:fB+5),fB-2);
% NBW = 1.5/N;
% Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / 3;
% plot(log10(f),dbp(Sqq*NBW),'m','Linewidth',2); grid on;
% text(0.001, -90, sprintf('NBW = %4.1E x f_s ',NBW),'Hor','right');
% legend(sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR));

% subplot(2,1,2);
spec_v = fft(v_a.*ds_hann(N))/(N/4);
% plot(log10(f),dbv(spec_v(1:N/2+1))); hold on;
snr_v = calculateSNR(spec_v(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
% plot(log10(f),dbv(spec(1:N/2+1))); hold on;
% plot(log10(f),dbv(spec_lsli(1:N/2+1))); grid on;
snr_lsli(i-1) = calculateSNR(spec_lsli(3:fB+5),fB-2);

% legend(sprintf('2-bit Quantizer SNR = %4.1fdB @ OSR = %d',snr_v,OSR),sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),sprintf('Leslie Singh SNR = %4.1fdB @ OSR = %d',snr_lsli,OSR));
end
figure;
plot(2:20,snr_lsli-snr); grid on;

% N = 32768; fs = 1e6; f = (33*fs)/N;
% order = 2; OSR = N/66; nlev = 2;
% u = sin(2*pi*f.*(0:1/fs:N/fs));
% sdin = csvread('Q.csv',1,0);
% sdin(:,2) = 1.8.*(double(sdin(:,2) > 0.9));
% v_tmp = sdin(1:N,2)';
% v_tmp = downsample(sdin(:,2),287)';
% v_tmp = v_tmp(10:32778);
% yin1 = csvread('von.csv',1,0);
% yin2 = csvread('vip1.csv',1,0);
% y = yin1(1:N,2)';
% y = downsample(yin1(:,2),287)';
% y = y(10:32778);
% fx = -0.5*fs:fs/N:0.5*fs; fint = (0.5*N+floor(N*0.8*f/fs)):((0.5*N+1)+floor(N*1.2*f/fs)); 
% hwv = fftshift(abs(fft(v_tmp'.*hanning(N+1))))./(N+1);
% hwu = fftshift(abs(fft(u'.*hanning(N+1))))./(N+1);
% plot(log10(fx),20.*log10(hwu)); hold on
% plot(log10(fx),20.*log10(hwv)); grid on
end