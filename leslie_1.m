N = 32768; fs = 1e6; f = (33*fs)/N;
u = 3*sin(2*pi*f.*(0:1/fs:N/fs)); w = (10e-4).*normrnd(0,1,1,N+1); 
order = 2; OSR = N/66; nlev = 2;
ntf = synthesizeNTF(order,OSR); plotPZ(ntf);
fx = -0.5*fs:fs/N:0.5*fs; fint = (0.5*N+3):((0.5*N+1)+ N*1.5*f/fs);

[snr,v_tmp,y] = sd_test_func(N,fs,f,0.7721.*u./3,w,order,OSR,ntf,nlev);
v_tmp12 = zeros(numel(v_tmp),12);
for i=1:numel(v_tmp)
    if (v_tmp(i) == -1)
        v_tmp12(i,:) = [1 zeros(1,11)];
    end
    if (v_tmp(i) == 1)
        v_tmp12(i,:) = [0 ones(1,11)];
    end
end
v = ds_quantize(3.*y,4); v_12 = zeros(numel(v),12);
for i=1:numel(v)
    if (v(i) == -3)
        v_12(i,:) = [1 zeros(1,11)];
    end
    if (v(i) == -1)
        v_12(i,:) = [1 1 zeros(1,10)];
    end
    if (v(i) == 1)
        v_12(i,:) = [0 1 zeros(1,10)];
    end
    if (v(i) == 3)
        v_12(i,:) = [0 ones(1,11)];
    end
end
v_tmp_1 = zeros(numel(v),12); v_tmp_2 = zeros(numel(v),12); v_tmp_3 = zeros(numel(v),12);
v_1 = zeros(numel(v),12); v_2 = zeros(numel(v),12); v_lsli = zeros(numel(v),12);
v_tmp_1(1,:) = bin_add(v_tmp12(1,:),zeros(1,12));
v_1(1,:) = bin_add(v_12(1,:),zeros(1,12));
v_tmp_2(1,:) = bin_add(v_tmp_1(1,:),zeros(1,12));
v_2(1,:) = bin_add(v_1(1,:),zeros(1,12));
v_tmp_3(1,:) = bin_add(sgn_xtnd(v_tmp12(1,:)),two_cmpl(sgn_xtnd(v_tmp_2(1,:))));
v_lsli(1,:) = bin_add(sgn_xtnd(v_2(1,:)),v_tmp_3(1,:));
for i=2:numel(v)
v_tmp_1(i,:) = bin_add(sgn_xtnd(v_tmp12(i,:)),two_cmpl(sgn_xtnd(v_tmp12(i-1,:))));
v_1(i,:) = bin_add(sgn_xtnd(v_12(i,:)),two_cmpl(sgn_xtnd(v_12(i-1,:))));
v_tmp_2(i,:) = bin_add(sgn_xtnd(v_tmp_1(i,:)),two_cmpl(sgn_xtnd(v_tmp_1(i-1,:))));
v_2(i,:) = bin_add(sgn_xtnd(v_1(i,:)),two_cmpl(sgn_xtnd(v_1(i-1,:))));
v_tmp_3(i,:) = bin_add(sgn_xtnd(v_tmp12(i,:)),two_cmpl(sgn_xtnd(v_tmp_2(i,:))));
v_lsli(i,:) = bin_add(sgn_xtnd(v_2(i,:)),v_tmp_3(i,:));
end
v_lsli_a = zeros(1,numel(v));
for i=1:numel(v)
v_lsli_a(i) = -v_lsli(i,1) + bi2de([0 v_lsli(i,2:12)],'left-msb')/(2^11);
end
v_lsli_a = 2.*v_lsli_a;
% bi2de(de2bi(floor(bi2de(v_tmp_3,'left-msb')*0.2165),12,'left-msb'),'left-msb')/(2^11)

[vw_tmp,xnw,xmaxw,yw] = simulateDSM(w,ntf,nlev);
vw_tmp12 = zeros(numel(vw_tmp),12);
for i=1:numel(vw_tmp)
    if (vw_tmp(i) == -1)
        vw_tmp12(i,:) = [1 zeros(1,11)];
    end
    if (vw_tmp(i) == 1)
        vw_tmp12(i,:) = [0 ones(1,11)];
    end
end
vw = ds_quantize(3.*yw,4); vw_12 = zeros(numel(vw),12);
for i=1:numel(vw)
    if (vw(i) == -2)
        vw_12(i,:) = [1 zeros(1,11)];
    end
    if (vw(i) == -1)
        vw_12(i,:) = [1 1 zeros(1,10)];
    end
    if (vw(i) == 1)
        vw_12(i,:) = [0 1 zeros(1,10)];
    end
    if (vw(i) == 2)
        vw_12(i,:) = [0 ones(1,11)];
    end
end
vw_tmp_1 = zeros(numel(vw),12); vw_tmp_2 = zeros(numel(vw),12); vw_tmp_3 = zeros(numel(vw),12);
vw_1 = zeros(numel(vw),12); vw_2 = zeros(numel(vw),12); vw_lsli = zeros(numel(vw),12);
vw_tmp_1(1,:) = bin_add(vw_tmp12(1,:),zeros(1,12));
vw_1(1,:) = bin_add(vw_12(1,:),zeros(1,12));
vw_tmp_2(1,:) = bin_add(vw_tmp_1(1,:),zeros(1,12));
vw_2(1,:) = bin_add(vw_1(1,:),zeros(1,12));
vw_tmp_3(1,:) = bin_add(sgn_xtnd(vw_tmp12(1,:)),two_cmpl(sgn_xtnd(vw_tmp_2(1,:))));
vw_lsli(1,:) = bin_add(sgn_xtnd(vw_2(1,:)),vw_tmp_3(1,:));
for i=2:numel(vw)
vw_tmp_1(i,:) = bin_add(sgn_xtnd(vw_tmp12(i,:)),two_cmpl(sgn_xtnd(vw_tmp12(i-1,:))));
vw_1(i,:) = bin_add(sgn_xtnd(vw_12(i,:)),two_cmpl(sgn_xtnd(vw_12(i-1,:))));
vw_tmp_2(i,:) = bin_add(sgn_xtnd(vw_tmp_1(i,:)),two_cmpl(sgn_xtnd(vw_tmp_1(i-1,:))));
vw_2(i,:) = bin_add(sgn_xtnd(vw_1(i,:)),two_cmpl(sgn_xtnd(vw_1(i-1,:))));
vw_tmp_3(i,:) = bin_add(sgn_xtnd(vw_tmp12(i,:)),two_cmpl(sgn_xtnd(vw_tmp_2(i,:))));
vw_lsli(i,:) = bin_add(sgn_xtnd(vw_2(i,:)),vw_tmp_3(i,:));
end
vw_lsli_a = zeros(1,numel(vw));
for i=1:numel(v)
vw_lsli_a(i) = -vw_lsli(i,1) + bi2de([0 vw_lsli(i,2:12)],'left-msb')/(2^11);
end
vw_lsli_a = 2.*vw_lsli_a;

hwv = fftshift(abs(fft(v'.*hanning(N+1))))./(N+1); hwvw = fftshift(abs(fft(vw'.*hanning(N+1))))./(N+1);
hwv_tmp = fftshift(abs(fft(v_tmp'.*hanning(N+1))))./(N+1); hwvw_tmp = fftshift(abs(fft(vw_tmp'.*hanning(N+1))))./(N+1);
hwv_lsli_a = fftshift(abs(fft(v_lsli_a'.*hanning(N+1))))./(N+1); hwvw_lsli_a = fftshift(abs(fft(vw_lsli_a'.*hanning(N+1))))./(N+1);

figure;
snrv_tmp = 10*log10(sum(hwv_tmp(fint).^2)/sum(hwvw_tmp(fint).^2));
snrv_lsli_a = 10*log10(sum(hwv_lsli_a(fint).^2)/sum(hwvw_lsli_a(fint).^2));
plot(log10(fx),20.*log10(hwv_tmp)); hold on; grid on;
snrv = 10*log10(sum(hwv(fint).^2)/sum(hwvw(fint).^2));
plot(log10(fx),20.*log10(hwv)); hold on; plot(log10(fx),20.*log10(hwv_lsli_a));
legend(sprintf('snrvtmp = %f',snrv_tmp),sprintf('snrv = %f',snrv),sprintf('snrvlslia = %f',snrv_lsli_a));


