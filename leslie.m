N = 32768; fs = 1e6; f = (33*fs)/N;
u = 3*sin(2*pi*f.*(0:1/fs:N/fs)); w = (10e-4).*normrnd(0,1,1,N+1); 
order = 2; OSR = N/66; nlev = 2;
ntf = synthesizeNTF(order,OSR); plotPZ(ntf);
fx = -0.5*fs:fs/N:0.5*fs; fint = (0.5*N+3):((0.5*N+1)+ N*1.5*f/fs);

[snr,v_tmp,y] = sd_test_func(N,fs,f,0.7721.*u./3,w,order,OSR,ntf,nlev);
v_tmp_1 = -0.5585.*[-3 -3 v_tmp(1:numel(v_tmp)-2)] + 0.775.*[-3 v_tmp(1:numel(v_tmp)-1)];
v = ds_quantize(3.*y,4);
v1a = v - [-3 v(1:numel(v)-1)]; v1a = v1a./4;
v2a = v1a - [-3 v1a(1:numel(v)-1)]; v2a = v2a./4;
v_lsli = v_tmp_1 + v2a;
% v = 0.5 + 0.5.*v;
% v1 = xor([1 1-v(1:numel(v)-1)],v); 
% v2 = xor([1 1-v1(1:numel(v1)-1)],v1); 
% v1 = 2.*v1 - 1; v2 = 2.*v2 - 1;

[vw_tmp,xnw,xmaxw,yw] = simulateDSM(w,ntf,nlev);
vw = ds_quantize(3.*yw,4);
% v1w = xor([1 1-vw(1:numel(vw)-1)],vw);
% v2w = xor([1 1-v1w(1:numel(v1w)-1)],v1w);
% v1w = 2.*v1w - 1; v2w = 2.*v2w - 1;
vw_tmp_1 = -0.5585.*[-3 -3 vw_tmp(1:numel(vw_tmp)-2)] + 0.775.*[-3 vw_tmp(1:numel(vw_tmp)-1)];
v1aw = vw - [-3 vw(1:numel(vw)-1)]; v1aw = v1aw./4;
v2aw = v1aw - [-3 v1aw(1:numel(v1aw)-1)]; v2aw = v2aw./4;
vw_lsli = vw_tmp_1 + v2aw;

hwv = fftshift(abs(fft(v'.*hanning(N+1))))./(N+1); hwvw = fftshift(abs(fft(vw'.*hanning(N+1))))./(N+1);
hwv_tmp = fftshift(abs(fft(v_tmp'.*hanning(N+1))))./(N+1); hwvw_tmp = fftshift(abs(fft(vw_tmp'.*hanning(N+1))))./(N+1);
hwv_tmp_1 = fftshift(abs(fft(v_tmp_1'.*hanning(N+1))))./(N+1); hwvw_tmp_1 = fftshift(abs(fft(vw_tmp_1'.*hanning(N+1))))./(N+1);
% hwv1 = fftshift(abs(fft(v1'.*hanning(N+1))))./(N+1); hwv1w = fftshift(abs(fft(v1w'.*hanning(N+1))))./(N+1);
% hwv2 = fftshift(abs(fft(v2'.*hanning(N+1))))./(N+1); hwv2w = fftshift(abs(fft(v2w'.*hanning(N+1))))./(N+1);
hwv1a = fftshift(abs(fft(v1a'.*hanning(N+1))))./(N+1); hwv1aw = fftshift(abs(fft(v1aw'.*hanning(N+1))))./(N+1);
hwv2a = fftshift(abs(fft(v2a'.*hanning(N+1))))./(N+1); hwv2aw = fftshift(abs(fft(v2aw'.*hanning(N+1))))./(N+1);
hwv_lsli = fftshift(abs(fft(v_lsli'.*hanning(N+1))))./(N+1); hwvw_lsli = fftshift(abs(fft(vw_lsli'.*hanning(N+1))))./(N+1);

figure;
snrv_tmp = 10*log10(sum(hwv_tmp(fint).^2)/sum(hwvw_tmp(fint).^2));
plot(log10(fx),20.*log10(hwv_tmp)); hold on;
snrv = 10*log10(sum(hwv(fint).^2)/sum(hwvw(fint).^2));
plot(log10(fx),20.*log10(hwv)); hold on; 
% snrv1 = 10*log10(sum(hwv1(fint).^2)/sum(hwv1w(fint).^2));
% plot(log10(fx),20.*log10(hwv1)); hold on; 
% snrv2 = 10*log10(sum(hwv2(fint).^2)/sum(hwv2w(fint).^2));
% plot(log10(fx),20.*log10(hwv2)); hold on; 
snrv1a = 10*log10(sum(hwv1a(fint).^2)/sum(hwv1aw(fint).^2));
plot(log10(fx),20.*log10(hwv1a)); hold on; 
snrv2a = 10*log10(sum(hwv2a(fint).^2)/sum(hwv2aw(fint).^2));
plot(log10(fx),20.*log10(hwv2a)); grid on;
legend(sprintf('snrvtmp = %f',snrv_tmp),sprintf('snrv = %f',snrv),sprintf('snrv1a = %f',snrv1a),sprintf('snrv2a = %f',snrv2a));

figure;
plot(log10(fx),20.*log10(hwv_tmp)); hold on; 
% snrv1 = 10*log10(sum(hwv1(fint).^2)/sum(hwv1w(fint).^2));
% plot(log10(fx),20.*log10(hwv1)); hold on; 
% snrv2 = 10*log10(sum(hwv2(fint).^2)/sum(hwv2w(fint).^2));
% plot(log10(fx),20.*log10(hwv2)); hold on; 
snrv_tmp_1 = 10*log10(sum(hwv_tmp_1(fint).^2)/sum(hwvw_tmp_1(fint).^2));
plot(log10(fx),20.*log10(hwv_tmp_1)); grid on; 
legend(sprintf('snrvtmp = %f',snrv_tmp),sprintf('snrvtmp_1 = %f',snrv_tmp_1));

figure;
snrv_lsli = 10*log10(sum(hwv_lsli(fint).^2)/sum(hwvw_lsli(fint).^2));
plot(log10(fx),20.*log10(hwv_tmp)); hold on; plot(log10(fx),20.*log10(hwv_lsli)); grid on;
legend(sprintf('snrvtmp = %f',snrv_tmp),sprintf('snrvlsli = %f',snrv_lsli));