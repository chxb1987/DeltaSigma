N = 32768; fs = 1e6; f = (33*fs)/N;
u = 3*sin(2*pi*f.*(0:1/fs:N/fs)); w = (10e-4).*normrnd(0,1,1,N+1);
order = 2; OSR = N/66; nlev = 2;
ntf = synthesizeNTF(order,OSR); plotPZ(ntf);
fx = -0.5*fs:fs/N:0.5*fs; fint = (0.5*N+3):((0.5*N+1)+ N*1.5*f/fs);

[snr,v_tmp,y,yt,vt_tmp] = sd_test_func(N,fs,f,0.7071.*u./3,w,order,OSR,ntf,nlev);
v_tmp16 = zeros(numel(v_tmp),16);
for i=1:numel(v_tmp)
    if (v_tmp(i) == -1)
%         v_tmp16(i,:) = [1 zeros(1,11)];
          v_tmp16(i,:) = [1 0 repmat([1 0],1,7)];
    end
    if (v_tmp(i) == 1)
%         v_tmp16(i,:) = [0 ones(1,11)];
          v_tmp16(i,:) = [0 1 repmat([0 1],1,7)];
    end
end
v = ds_quantize(3.*y,4); v_16 = zeros(numel(v),16);
for i=1:numel(v)
    if (v(i) == -3)
        v_16(i,:) = [1 zeros(1,15)];
    end
    if (v(i) == -1)
%         v_16(i,:) = [1 1 zeros(1,10)];
          v_16(i,:) = [1 1 repmat([0 1],1,7)];
    end
    if (v(i) == 1)
%         v_16(i,:) = [0 1 zeros(1,10)];
          v_16(i,:) = [0 0 repmat([1 0],1,7)];
    end
    if (v(i) == 3)
        v_16(i,:) = [0 ones(1,15)];
    end
end
v_tmp_1 = zeros(numel(v),16); v_tmp_2 = zeros(numel(v),16); v_tmp_3 = zeros(numel(v),16);
v_1 = zeros(numel(v),16); v_2 = zeros(numel(v),16); v_lsli = zeros(numel(v),16);
v_tmp_1(1,:) = bin_add(v_tmp16(1,:),zeros(1,16));
v_1(1,:) = bin_add(v_16(1,:),zeros(1,16));
v_tmp_2(1,:) = bin_add(v_tmp_1(1,:),zeros(1,16));
v_2(1,:) = bin_add(v_1(1,:),zeros(1,16));
v_tmp_3(1,:) = bin_add(sgn_xtnd(v_tmp16(1,:)),two_cmpl(sgn_xtnd(v_tmp_2(1,:))));
v_lsli(1,:) = bin_add(v_2(1,:),v_tmp_3(1,:));
for i=2:numel(v)
v_tmp_1(i,:) = bin_add(sgn_xtnd(v_tmp16(i,:)),two_cmpl(sgn_xtnd(v_tmp16(i-1,:))));
v_1(i,:) = bin_add(sgn_xtnd(v_16(i,:)),two_cmpl(sgn_xtnd(v_16(i-1,:))));
v_tmp_2(i,:) = bin_add(sgn_xtnd(v_tmp_1(i,:)),two_cmpl(sgn_xtnd(v_tmp_1(i-1,:))));
v_2(i,:) = bin_add(sgn_xtnd(v_1(i,:)),two_cmpl(sgn_xtnd(v_1(i-1,:))));
v_tmp_3(i,:) = bin_add(sgn_xtnd(v_tmp16(i,:)),two_cmpl(sgn_xtnd(v_tmp_2(i,:))));
v_lsli(i,:) = bin_add(v_2(i,:),v_tmp_3(i,:));
end
v_lsli_a = zeros(1,numel(v));
for i=1:numel(v)
v_lsli_a(i) = -v_lsli(i,1) + bi2de([0 v_lsli(i,2:16)],'left-msb')/(2^15);
end
v_lsli_a = 2.*v_lsli_a;
% bi2de(de2bi(floor(bi2de(v_tmp_3,'left-msb')*0.2165),16,'left-msb'),'left-msb')/(2^15)

[vw_tmp,xnw,xmaxw,yw] = simulateDSM(w,ntf,nlev);
vw_tmp16 = zeros(numel(vw_tmp),16);
for i=1:numel(vw_tmp)
    if (vw_tmp(i) == -1)
%         vw_tmp16(i,:) = [1 zeros(1,11)];
          v_tmp16(i,:) = [1 0 repmat([1 0],1,7)];
    end
    if (vw_tmp(i) == 1)
%         vw_tmp16(i,:) = [0 ones(1,11)];
          v_tmp16(i,:) = [0 1 repmat([0 1],1,7)];
    end
end
vw = ds_quantize(3.*yw,4); vw_16 = zeros(numel(vw),16);
for i=1:numel(vw)
    if (vw(i) == -3)
        vw_16(i,:) = [1 zeros(1,15)];
    end
    if (vw(i) == -1)
%         vw_16(i,:) = [1 1 zeros(1,10)];
          v_16(i,:) = [1 1 repmat([0 1],1,7)];
    end
    if (vw(i) == 1)
%         vw_16(i,:) = [0 1 zeros(1,10)];
          v_16(i,:) = [0 0 repmat([1 0],1,7)];
    end
    if (vw(i) == 3)
        vw_16(i,:) = [0 ones(1,15)];
    end
end
vw_tmp_1 = zeros(numel(vw),16); vw_tmp_2 = zeros(numel(vw),16); vw_tmp_3 = zeros(numel(vw),16);
vw_1 = zeros(numel(vw),16); vw_2 = zeros(numel(vw),16); vw_lsli = zeros(numel(vw),16);
vw_tmp_1(1,:) = bin_add(vw_tmp16(1,:),zeros(1,16));
vw_1(1,:) = bin_add(vw_16(1,:),zeros(1,16));
vw_tmp_2(1,:) = bin_add(vw_tmp_1(1,:),zeros(1,16));
vw_2(1,:) = bin_add(vw_1(1,:),zeros(1,16));
vw_tmp_3(1,:) = bin_add(sgn_xtnd(vw_tmp16(1,:)),two_cmpl(sgn_xtnd(vw_tmp_2(1,:))));
vw_lsli(1,:) = bin_add(vw_2(1,:),vw_tmp_3(1,:));
for i=2:numel(vw)
vw_tmp_1(i,:) = bin_add(sgn_xtnd(vw_tmp16(i,:)),two_cmpl(sgn_xtnd(vw_tmp16(i-1,:))));
vw_1(i,:) = bin_add(sgn_xtnd(vw_16(i,:)),two_cmpl(sgn_xtnd(vw_16(i-1,:))));
vw_tmp_2(i,:) = bin_add(sgn_xtnd(vw_tmp_1(i,:)),two_cmpl(sgn_xtnd(vw_tmp_1(i-1,:))));
vw_2(i,:) = bin_add(sgn_xtnd(vw_1(i,:)),two_cmpl(sgn_xtnd(vw_1(i-1,:))));
vw_tmp_3(i,:) = bin_add(sgn_xtnd(vw_tmp16(i,:)),two_cmpl(sgn_xtnd(vw_tmp_2(i,:))));
vw_lsli(i,:) = bin_add(vw_2(i,:),vw_tmp_3(i,:));
end
vw_lsli_a = zeros(1,numel(vw));
for i=1:numel(v)
vw_lsli_a(i) = -vw_lsli(i,1) + bi2de([0 vw_lsli(i,2:16)],'left-msb')/(2^15);
end
vw_lsli_a = 2.*vw_lsli_a;

vt_tmp16 = zeros(numel(vt_tmp),16);
for i=1:numel(vt_tmp)
    if (vt_tmp(i) == -1)
%         vt_tmp16(i,:) = [1 zeros(1,11)];
          vt_tmp16(i,:) = [1 0 repmat([1 0],1,7)];
    end
    if (vt_tmp(i) == 1)
%         vt_tmp16(i,:) = [0 ones(1,11)];
          vt_tmp16(i,:) = [0 1 repmat([0 1],1,7)];
    end
end
vt = ds_quantize(3.*yt,4); vt_16 = zeros(numel(vt),16);
for i=1:numel(vt)
    if (vt(i) == -3)
        vt_16(i,:) = [1 zeros(1,15)];
    end
    if (vt(i) == -1)
%         vt_16(i,:) = [1 1 zeros(1,10)];
          vt_16(i,:) = [1 1 repmat([0 1],1,7)];
    end
    if (vt(i) == 1)
%         vt_16(i,:) = [0 1 zeros(1,10)];
          vt_16(i,:) = [0 0 repmat([1 0],1,7)];
    end
    if (vt(i) == 3)
        vt_16(i,:) = [0 ones(1,15)];
    end
end
vt_tmp_1 = zeros(numel(vt),16); vt_tmp_2 = zeros(numel(vt),16); vt_tmp_3 = zeros(numel(vt),16);
vt_1 = zeros(numel(vt),16); vt_2 = zeros(numel(vt),16); vt_lsli = zeros(numel(vt),16);
vt_tmp_1(1,:) = bin_add(vt_tmp16(1,:),zeros(1,16));
vt_1(1,:) = bin_add(vt_16(1,:),zeros(1,16));
vt_tmp_2(1,:) = bin_add(vt_tmp_1(1,:),zeros(1,16));
vt_2(1,:) = bin_add(vt_1(1,:),zeros(1,16));
vt_tmp_3(1,:) = bin_add(sgn_xtnd(vt_tmp16(1,:)),two_cmpl(sgn_xtnd(vt_tmp_2(1,:))));
vt_lsli(1,:) = bin_add(vt_2(1,:),vt_tmp_3(1,:));
for i=2:numel(vt)
vt_tmp_1(i,:) = bin_add(sgn_xtnd(vt_tmp16(i,:)),two_cmpl(sgn_xtnd(vt_tmp16(i-1,:))));
vt_1(i,:) = bin_add(sgn_xtnd(vt_16(i,:)),two_cmpl(sgn_xtnd(vt_16(i-1,:))));
vt_tmp_2(i,:) = bin_add(sgn_xtnd(vt_tmp_1(i,:)),two_cmpl(sgn_xtnd(vt_tmp_1(i-1,:))));
vt_2(i,:) = bin_add(sgn_xtnd(vt_1(i,:)),two_cmpl(sgn_xtnd(vt_1(i-1,:))));
vt_tmp_3(i,:) = bin_add(sgn_xtnd(vt_tmp16(i,:)),two_cmpl(sgn_xtnd(vt_tmp_2(i,:))));
vt_lsli(i,:) = bin_add(vt_2(i,:),vt_tmp_3(i,:));
end
vt_lsli_a = zeros(1,numel(vt));
for i=1:numel(vt)
vt_lsli_a(i) = -vt_lsli(i,1) + bi2de([0 vt_lsli(i,2:16)],'left-msb')/(2^15);
end
vt_lsli_a = 2.*vt_lsli_a;

hwv = fftshift(abs(fft(v'.*hanning(N+1))))./(N+1); hwvw = fftshift(abs(fft(vw'.*hanning(N+1))))./(N+1);
hwv_tmp = fftshift(abs(fft(v_tmp'.*hanning(N+1))))./(N+1); hwvw_tmp = fftshift(abs(fft(vw_tmp'.*hanning(N+1))))./(N+1);
hwv_lsli_a = fftshift(abs(fft(v_lsli_a'.*hanning(N+1))))./(N+1); hwvw_lsli_a = fftshift(abs(fft(vw_lsli_a'.*hanning(N+1))))./(N+1);

hwvt = fftshift(abs(fft(vt'.*hanning(N+1))))./(N+1); 
hwvt_tmp = fftshift(abs(fft(vt_tmp'.*hanning(N+1))))./(N+1); 
hwvt_lsli_a = fftshift(abs(fft(vt_lsli_a'.*hanning(N+1))))./(N+1);

figure;
snrv_tmp = 10*log10(sum(hwv_tmp(fint).^2)/sum(hwvw_tmp(fint).^2));
snrv_lsli_a = 10*log10(sum(hwv_lsli_a(fint).^2)/sum(hwvw_lsli_a(fint).^2));
plot(log10(fx),20.*log10(hwvw_tmp)); hold on; grid on;
snrv = 10*log10(sum(hwv(fint).^2)/sum(hwvw(fint).^2));
plot(log10(fx),20.*log10(hwvw)); hold on; plot(log10(fx),20.*log10(hwvw_lsli_a));
legend(sprintf('snrvtmp = %f',snrv_tmp),sprintf('snrv = %f',snrv),sprintf('snrvlslia = %f',snrv_lsli_a));


