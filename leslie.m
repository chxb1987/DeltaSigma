N = 32768; OSR = N/66; nlev = 2;
H = synthesizeNTF(2,OSR,1);
fB = ceil(N/(2*OSR)); ftest = floor(2/3*fB);
u = 0.5*sin(2*pi*ftest/N*(0:N-1));	% half-scale sine-wave input
[v_tmp,xn,xnmax,y] = simulateDSM(u,H,nlev); 
f = linspace(0,0.5,N/2+1);

% figure;
% plot(v_tmp); hold on; plot(y); plot(2.*u); grid on;
v_tmp12 = zeros(numel(v_tmp),12);
for i=1:numel(v_tmp)
    if (v_tmp(i) == -1)
        v_tmp12(i,:) = [1 1 zeros(1,10)];
%         v_tmp12(i,:) = [1 zeros(1,11)];
%         v_tmp12(i,:) = [1 0 repmat([1 0],1,5)];
    end
    if (v_tmp(i) == 1)
        v_tmp12(i,:) = [0 1 zeros(1,10)];
%         v_tmp12(i,:) = [0 ones(1,11)];
%         v_tmp12(i,:) = [0 1 repmat([0 1],1,5)];
    end
end

v_tmp_a = zeros(1,numel(v_tmp));
for i=1:numel(v_tmp)
v_tmp_a(i) = -v_tmp12(i,1) + bi2de([0 v_tmp12(i,2:12)],'left-msb')/(2^11);
end
% v_tmp_a = 2.*v_tmp_a;
% figure;
% subplot(2,1,1); plot(v_tmp); hold on; plot(y); plot(2.*u); grid on;
% subplot(2,1,2); plot(v_tmp); hold on; plot(v_tmp_a); grid on;

v = ds_quantize(7.*y,8); v_12 = zeros(numel(v),12);
for i=1:numel(v)
    if (v(i) == -7)
        v_12(i,:) = [1 0 0 1 zeros(1,8)];
%         v_12(i,:) = [1 zeros(1,11)];
    end
    if (v(i) == -5)
        v_12(i,:) = [1 0 1 1 zeros(1,8)];
%         v_12(i,:) = [1 1 zeros(1,10)];
%         v_12(i,:) = [1 1 repmat([0 1],1,5)];
    end
    if (v(i) == -3)
        v_12(i,:) = [1 1 0 1 zeros(1,8)];
%         v_12(i,:) = [1 zeros(1,11)];
    end
    if (v(i) == -1)
        v_12(i,:) = [1 1 1 1 zeros(1,8)];
%         v_12(i,:) = [1 1 zeros(1,10)];
%         v_12(i,:) = [1 1 repmat([0 1],1,5)];
    end
    if (v(i) == 1)
        v_12(i,:) = [0 0 0 1 zeros(1,8)];
%         v_12(i,:) = [0 1 zeros(1,10)];
%         v_12(i,:) = [0 0 repmat([1 0],1,5)];
    end
    if (v(i) == 3)
        v_12(i,:) = [0 0 1 1 zeros(1,8)];
%         v_12(i,:) = [0 ones(1,11)];
    end
    if (v(i) == 5)
        v_12(i,:) = [0 1 0 1 zeros(1,8)];
%         v_12(i,:) = [1 zeros(1,11)];
    end
    if (v(i) == 7)
        v_12(i,:) = [0 1 1 1 zeros(1,8)];
%         v_12(i,:) = [1 1 zeros(1,10)];
%         v_12(i,:) = [1 1 repmat([0 1],1,5)];
    end
end

v_a = zeros(1,numel(v_tmp));
for i=1:numel(v_tmp)
v_a(i) = -v_12(i,1) + bi2de([0 v_12(i,2:12)],'left-msb')/(2^11);
end
% v_a = (8/7).*v_a;
% figure;
% subplot(2,1,1); plot(v_tmp); hold on; plot(y); plot(2.*u); grid on;
% subplot(2,1,2); plot(v_a); hold on; plot(6.*u); grid on;

v_tmp_1 = zeros(numel(v),12); v_1 = zeros(numel(v),12); v_2 = zeros(numel(v),12); v_lsli = zeros(numel(v),12);

v_tmp_1(1,:) = zeros(1,12);
v_1(1,:) = bin_add(v_12(1,:),zeros(1,12));
v_2(1,:) = bin_add(v_1(1,:),zeros(1,12));
v_lsli(1,:) = bin_add(v_2(1,:),v_tmp_1(1,:));

v_tmp_1(2,:) = bin_add(sgn_xtnd(sgn_xtnd(v_tmp12(1,:))),sgn_xtnd(v_tmp12(1,:)));
v_1(2,:) = bin_add(sgn_xtnd(v_12(2,:)),two_cmpl(sgn_xtnd(v_12(1,:))));
v_2(2,:) = bin_add(sgn_xtnd(v_1(2,:)),two_cmpl(sgn_xtnd(v_1(1,:))));
v_lsli(2,:) = bin_add(v_2(2,:),v_tmp_1(2,:));

for i=3:numel(v)
v_tmp_1(i,:) = bin_add(bin_add(sgn_xtnd(sgn_xtnd(v_tmp12(i-1,:))),sgn_xtnd(v_tmp12(i-1,:))),two_cmpl(sgn_xtnd(v_tmp12(i-2,:))));
v_1(i,:) = bin_add(sgn_xtnd(v_12(i,:)),two_cmpl(sgn_xtnd(v_12(i-1,:))));
v_2(i,:) = bin_add(sgn_xtnd(v_1(i,:)),two_cmpl(sgn_xtnd(v_1(i-1,:))));
v_lsli(i,:) = bin_add(v_2(i,:),v_tmp_1(i,:));
end

v_tmp_1_a = zeros(1,numel(v));
for i=1:numel(v)
v_tmp_1_a(i) = -v_tmp_1(i,1) + bi2de([0 v_tmp_1(i,2:12)],'left-msb')/(2^11);
end
v_1_a = zeros(1,numel(v));
for i=1:numel(v)
v_1_a(i) = -v_1(i,1) + bi2de([0 v_1(i,2:12)],'left-msb')/(2^11);
end
v_2_a = zeros(1,numel(v));
for i=1:numel(v)
v_2_a(i) = -v_2(i,1) + bi2de([0 v_2(i,2:12)],'left-msb')/(2^11);
end
v_lsli_a = zeros(1,numel(v));
for i=1:numel(v)
v_lsli_a(i) = -v_lsli(i,1) + bi2de([0 v_lsli(i,2:12)],'left-msb')/(2^11);
end
% v_lsli_a = (64).*v_lsli_a;
% figure;
% plot(v_lsli_a); hold on; plot(2.*u); grid on;

spec = fft(v_tmp_a.*ds_hann(N))/(N/4);
figure;
subplot(2,1,1);
plot(log10(f),dbv(spec(1:N/2+1))); hold on;
snr = calculateSNR(spec(3:fB+1),ftest-2);
NBW = 1.5/N;
Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / 3;
plot(log10(f),dbp(Sqq*NBW),'m','Linewidth',2); grid on;
text(0.001, -90, sprintf('NBW = %4.1E x f_s ',NBW),'Hor','right');
legend(sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),'1-bit Expected');

subplot(2,1,2);
spec_v = fft(v_a.*ds_hann(N))/(N/4);
plot(log10(f),dbv(spec_v(1:N/2+1))); hold on;
snr_v = calculateSNR(spec_v(3:fB+1),ftest-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
plot(log10(f),dbv(spec(1:N/2+1))); hold on;
plot(log10(f),dbv(spec_lsli(1:N/2+1))); grid on;
snr_lsli = calculateSNR(spec_lsli(3:fB+1),ftest-2);

legend(sprintf('2-bit Quantizer SNR = %4.1fdB @ OSR = %d',snr_v,OSR),sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),sprintf('Leslie Singh SNR = %4.1fdB @ OSR = %d',snr_lsli,OSR));
