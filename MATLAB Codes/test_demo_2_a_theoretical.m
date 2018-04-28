function [] = test_demo_2_a_theoretical(flag)
if (flag == 0)
N = 32768; OSR = N/66; nlev = 2; f = linspace(0,0.5,N/2+1);
z = zpk('z',1);
H = ((z-1)^2)/(z^2-1.567*z+0.7835);
fB = ceil(N/(2*OSR));

u_theo = 0.5*sin(2*pi*fB/N*(0:N-1));	
[vtheo,~,~,y_theo] = simulateDSM(u_theo,H,nlev);    
vtheo_a = vtheo./2;

pso_gfunc_lsli([1,1.567,0.7835],vtheo_a,y_theo,1)
idl = zeros(1,19);
for i=1:numel(idl)
if (i == 1)
idl(i) = 20*log10(2^(i+1)-1) - 20*log10(2^(i)-1);
else
idl(i) = idl(i-1) + 20*log10(2^(i+1)-1) - 20*log10(2^(i)-1);
end
end
hold on; plot(2:20,idl,'-.k');
title('SNR Improvement vs Number of Quantizer Bits for Leslie-Singh Topology')
xlabel('B (Quantizer Bits)')
ylabel('SNR Improvement (dB)')
legend('Improvement for Leslie-Singh Topology','Theoretical Improvement for M-bit, 2nd order Delta-Sigma Modulator')

optim = [1,1.567,0.7835];
v = ds_quantize(((2^(4))-1).*y_theo,2^(4)); 
v_a = v./(2^(4));
v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];
    
v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;
spec = fft(vtheo_a.*ds_hann(N))/(N/4);
snr = calculateSNR(spec(3:fB+5),fB-2);

specv = fft(v_a.*ds_hann(N))/(N/4);
snrv = calculateSNR(specv(3:fB+5),fB-2);
  
spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);
    
figure;
subplot(2,1,1);
plot(log10(f),dbv(spec(1:N/2+1)),'k'); grid on; legend(sprintf('1-bit Delta-Sigma Modulator Output SNR = %4.1fdB @ OSR = %d',snr,OSR)); 
subplot(2,1,2);
plot(log10(f),dbv(specv(1:N/2+1)),'k'); grid on;
legend(sprintf('4-bit Quantized Delta-Sigma Modulator Intermediate Signal SNR = %4.1fdB @ OSR = %d',snrv,OSR));
   
figure;
subplot(2,1,1);
plot(log10(f),dbv(spec(1:N/2+1)),'k'); grid on; legend(sprintf('1-bit Delta-Sigma Modulator Output SNR = %4.1fdB @ OSR = %d',snr,OSR));
subplot(2,1,2);
plot(log10(f),dbv(spec_lsli(1:N/2+1)),'k'); grid on;
legend(sprintf('4-bit Leslie-Singh Topology Output SNR = %4.1fdB @ OSR = %d',snr_lsli,OSR));
end

if (flag == 1)
N = 32768; OSR = N/66; nlev = 2; f = linspace(0,0.5,N/2+1);
z = zpk('z',1);
H = ((z-1)^2)/(z^2-1.567*z+0.7835);
fB = ceil(N/(2*OSR));
A = 0.1:0.001:1; interm = zeros(1,numel(A));

for i=1:numel(A)
u_theo = A(i)*sin(2*pi*fB/N*(0:N-1));	% half-scale sine-wave input
[vtheo,~,~,y_theo] = simulateDSM(u_theo,H,nlev); 
vtheo_a = vtheo./2;

optim = [1,1.567,0.7835];
v = ds_quantize(((2^(4))-1).*y_theo,2^(4)); 
v_a = v./(2^(4));
v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];
    
v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;
spec = fft(vtheo_a.*ds_hann(N))/(N/4);
snr = calculateSNR(spec(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);
interm(i) = snr_lsli - snr;
end
figure;
plot(A,interm,'k'); grid on;
title('SNR Improvement vs Input Signal Amplitude for 4-bit Leslie-Singh Topology')
xlabel('Amplitude')
ylabel('SNR Improvement (dB)')
end

if (flag == 2)
% N = 2.^(10:20); OSR = N./66; nlev = 2; 
m = 3:2:33; N = 32768.*ones(1,numel(m)); OSR = N(1)./(2.*m); nlev = 2; 
z = zpk('z',1);
H = ((z-1)^2)/(z^2-1.567*z+0.7835);
interm = zeros(1,numel(N));

for i=1:numel(N)
fB = ceil(N(i)/(2*OSR(i)));
u_theo = 0.6*sin(2*pi*fB/N(i)*(0:N(i)-1));	% half-scale sine-wave input
[vtheo,~,~,y_theo] = simulateDSM(u_theo,H,nlev); 
vtheo_a = vtheo./2;

optim = [1,1.567,0.7835];
v = ds_quantize(((2^(4))-1).*y_theo,2^(4)); 
v_a = v./(2^(4));
v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];
    
v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;
spec = fft(vtheo_a.*ds_hann(N(i)))/(N(i)/4);
snr = calculateSNR(spec(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N(i)))/(N(i)/4);
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);
interm(i) = snr_lsli - snr;
end
figure;
plot(OSR,interm,'k'); grid on;
title('SNR Improvement vs OSR for 4-bit Leslie-Singh Topology')
xlabel('OSR')
ylabel('SNR Improvement (dB)')
end

if (flag == 3)
N = 32768; OSR = N/66; nlev = 2; f = linspace(0,0.5,N/2+1);
z = zpk('z',1);
H = ((z-1)^2)/(z^2-1.567*z+0.7835);
fB = ceil(N/(2*OSR));

u_theo = 0.6*sin(2*pi*fB/N*(0:N-1));	% half-scale sine-wave input
[vtheo,~,~,y_theo] = simulateDSM(u_theo,H,nlev);    
vtheo_a = vtheo./2; 

% alpha = 0.4:0.01:0.6; beta = 0.95:0.005:1.05;
alpha = 0.95:0.00001:1.05; beta = 0.95:0.00001:1.05;
interma = zeros(1,numel(alpha)); intermb = zeros(1,numel(beta));

for i=1:numel(alpha)
optim = [1,2-0.433*alpha(i),1-0.433*alpha(i)+0.2165*alpha(i)^2];
v = ds_quantize(((2^(4))-1).*y_theo,2^(4)); 
v_a = v./(2^(4));
v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];
    
v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;
spec = fft(vtheo_a.*ds_hann(N))/(N/4);
snr = calculateSNR(spec(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);
interma(i) = snr_lsli - snr;
end
figure;
plot(alpha,interma,'k'); grid on;
title('SNR Improvement vs Variation in alpha for 4-bit Leslie-Singh Topology')
xlabel('Variation in alpha')
ylabel('SNR Improvement (dB)')

for i=1:numel(beta)
optim = [beta(i),2*beta(i)-0.433,beta(i)^2-0.433*beta(i)+0.2165];
v = ds_quantize(((2^(4))-1).*y_theo,2^(4)); 
v_a = v./(2^(4));
v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];
    
v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;
spec = fft(vtheo_a.*ds_hann(N))/(N/4);
snr = calculateSNR(spec(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);
intermb(i) = snr_lsli - snr;
end
figure;
plot(beta,intermb,'k'); grid on;
title('SNR Improvement vs Variation in beta for 4-bit Leslie-Singh Topology')
xlabel('Variation in beta')
ylabel('SNR Improvement (dB)')

da = diff(interma)./0.00001;
db = diff(intermb)./0.00001;
figure;
plot(alpha(2:end),da,'k'); grid on;
title('First Derivative of SNR Improvement vs Variation in alpha for 4-bit Leslie-Singh Topology')
xlabel('Variation in alpha')
ylabel('First Derivative of SNR Improvement (dB)')
figure;
plot(beta(2:end),db,'k'); grid on;
title('First Derivative of SNR Improvement vs Variation in beta for 4-bit Leslie-Singh Topology')
xlabel('Variation in beta')
ylabel('First Derivative of SNR Improvement (dB)')
end

if (flag == 4)
N = 32768; fs = 1e6;
gamma = 1e3;
theta = 6;
delta = 1;

b = [0 1];
a = [delta -delta];
[h,w] = freqz(b,a,N,'whole',fs);
figure;
plot(log10(w),20*log10(abs(h)),'k'); grid on; hold on;

alpha = gamma/(1+delta*(1+gamma));
beta = delta*(1+gamma)/(1+delta*(1+gamma));
b = [0 alpha];
a = [1 -beta];
[h,w] = freqz(b,a,N,'whole',fs);
plot(log10(w),20*log10(abs(h)),'--k'); 

b = (gamma*theta).*[0 1 1];
a = [((2*delta+2)*gamma+(1+delta)*theta+delta*gamma*theta) (theta-(2+4*delta)*gamma) delta*(2*gamma-(1+gamma)*theta)];
[h,w] = freqz(b,a,N,'whole',fs);
plot(log10(w),20*log10(abs(h)),'-.k');
title('Integrator Transfer Function Magnitude Response')
xlabel('log10(f)')
ylabel('20log10|H(z)|')
legend('Ideal','Non-ideal High First Pole Approximation, A(0) = 60 dB, UGF = 6 MHz','Non-ideal Comparable First Pole Approximation, A(0) = 60 dB, UGF = 6 MHz');
end

if (flag == 5)
N = 32768; OSR = N/66; nlev = 2; f = linspace(0,0.5,N/2+1);
z = zpk('z',1);
gamma = 1e3;
theta = 6;
delta = 1;
alpha = gamma/(1+delta*(1+gamma));
beta = delta*(1+gamma)/(1+delta*(1+gamma));
H = ((z-beta)^2)/(z^2-(2*beta-0.433*alpha)*z+(beta^2+0.2165*alpha^2-0.433*alpha*beta)); H1 = ((z-1)^2)/(z^2-1.567*z+0.7835);
fB = ceil(N/(2*OSR));
u_theo = 0.6*sin(2*pi*fB/N*(0:N-1));	% half-scale sine-wave input
[vtheo,~,~,y_theo] = simulateDSM(u_theo,H,nlev); [vtheo1,~,~,y_theo1] = simulateDSM(u_theo,H1,nlev);     
spec = fft(vtheo.*ds_hann(N))/(N/4); spec1 = fft(vtheo1.*ds_hann(N))/(N/4);
snr = calculateSNR(spec(3:fB+5),fB-2); snr1 = calculateSNR(spec1(3:fB+5),fB-2);
figure;
plot(log10(f),dbv(spec(1:N/2+1)),'-.k'); grid on; 
hold on;
plot(log10(f),dbv(spec1(1:N/2+1)),'k'); grid on; legend(sprintf('Non-ideal Integrator, SNR = %4.1fdB @ OSR = %d',snr,OSR),sprintf('Ideal Integrator, SNR = %4.1fdB @ OSR = %d',snr1,OSR)); 
end

if (flag == 6)
N = 32768; OSR = N/66; nlev = 2; f = linspace(0,0.5,N/2+1);
z = zpk('z',1);
H = ((z-1)^2)/(z^2-1.567*z+0.7835);
fB = ceil(N/(2*OSR));

u_theo = 0.6*sin(2*pi*fB/N*(0:N-1));	% half-scale sine-wave input
[vtheo,~,~,y_theo] = simulateDSM(u_theo,H,nlev);    
vtheo_a = vtheo./2; 

gm = 0:0.01:0.99; dt = 0:0.01:0.99;
intermgm = zeros(1,numel(gm)); intermdt = zeros(1,numel(dt));

for i=1:numel(gm)
optim = [1,2-0.433,1-0.433+0.433*gm(i)];
v = ds_quantize(((2^(4))-1).*y_theo,2^(4)); 
v_a = v./(2^(4));
v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];
    
v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;
spec = fft(vtheo_a.*ds_hann(N))/(N/4);
snr = calculateSNR(spec(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);
intermgm(i) = snr_lsli - snr;
end
figure;
plot(gm,intermgm,'k'); grid on; hold on
title('SNR Improvement vs Variation in gamma and delta for 4-bit Leslie-Singh Topology')
xlabel('Variation in gamma and delta')
ylabel('SNR Improvement (dB)')

for i=1:numel(dt)
optim = [1,2-dt(i),1-dt(i)+0.5*dt(i)];
v = ds_quantize(((2^(4))-1).*y_theo,2^(4)); 
v_a = v./(2^(4));
v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];
    
v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;
spec = fft(vtheo_a.*ds_hann(N))/(N/4);
snr = calculateSNR(spec(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);
intermdt(i) = snr_lsli - snr;
end
plot(dt,intermdt,'-.k'); 
legend('gamma','delta');
end
end
