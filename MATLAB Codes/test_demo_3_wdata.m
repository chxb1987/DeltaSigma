N = 32768; OSR = N/66; fB = ceil(N/(2*OSR)); 
sdin = csvread('Q_op3.csv',1,0);
sdin(:,2) = 2.*(double(sdin(:,2) > 0.9))-1;
v_tmp = sdin(10:(N+9),2)';

optim = [1,-57,-58];
z = zpk('z',1);
H = ((z-optim(1))^2)/(z^2-optim(2)*z+optim(3));

yin1 = csvread('von_op3.csv',1,0);
% yin2 = csvread('vip1_op3.csv',1,0);
yin1(:,2) = (yin1(:,2) - 0.9)./0.1;
% yin2(:,2) = (yin2(:,2) - 0.9)./0.125;
y = yin1(10:(N+9),2)';
f = linspace(0,0.5,N/2+1);

v_tmp16 = zeros(numel(v_tmp),16);
for i=1:numel(v_tmp)
    if (v_tmp(i) == -1)
        v_tmp16(i,:) = [ones(1,8) 1 zeros(1,7)];
%         v_tmp16(i,:) = [1 zeros(1,11)];
%         v_tmp16(i,:) = [1 0 repmat([1 0],1,5)];
    end
    if (v_tmp(i) == 1)
        v_tmp16(i,:) = [zeros(1,8) 1 zeros(1,7)];
%         v_tmp16(i,:) = [0 ones(1,11)];
%         v_tmp16(i,:) = [0 1 repmat([0 1],1,5)];
    end
end

v_tmp_a = zeros(1,numel(v_tmp));
for i=1:numel(v_tmp)
v_tmp_a(i) = -v_tmp16(i,1)*(2^7) + bi2de([0 v_tmp16(i,2:8)],'left-msb') + bi2de(v_tmp16(i,9:16),'left-msb')/(2^8);
end

v = ds_quantize(7.*y,8); v_16 = zeros(numel(v),16);
for i=1:numel(v)
    if (v(i) == -7)
        v_16(i,:) = [ones(1,8) 0 0 1 zeros(1,5)];
%         v_16(i,:) = [1 zeros(1,11)];
    end
    if (v(i) == -5)
        v_16(i,:) = [ones(1,8) 0 1 1 zeros(1,5)];
%         v_16(i,:) = [1 1 zeros(1,10)];
%         v_16(i,:) = [1 1 repmat([0 1],1,5)];
    end
    if (v(i) == -3)
        v_16(i,:) = [ones(1,8) 1 0 1 zeros(1,5)];
%         v_16(i,:) = [1 zeros(1,11)];
    end
    if (v(i) == -1)
        v_16(i,:) = [ones(1,8) 1 1 1 zeros(1,5)];
%         v_16(i,:) = [1 1 zeros(1,10)];
%         v_16(i,:) = [1 1 repmat([0 1],1,5)];
    end
    if (v(i) == 1)
        v_16(i,:) = [zeros(1,8) 0 0 1 zeros(1,5)];
%         v_16(i,:) = [0 1 zeros(1,10)];
%         v_16(i,:) = [0 0 repmat([1 0],1,5)];
    end
    if (v(i) == 3)
        v_16(i,:) = [zeros(1,8) 0 1 1 zeros(1,5)];
%         v_16(i,:) = [0 ones(1,11)];
    end
    if (v(i) == 5)
        v_16(i,:) = [zeros(1,8) 1 0 1 zeros(1,5)];
%         v_16(i,:) = [1 zeros(1,11)];
    end
    if (v(i) == 7)
        v_16(i,:) = [zeros(1,8) 1 1 1 zeros(1,5)];
%         v_16(i,:) = [1 1 zeros(1,10)];
%         v_16(i,:) = [1 1 repmat([0 1],1,5)];
    end
end

v_a = zeros(1,numel(v_tmp));
for i=1:numel(v_tmp)
v_a(i) = -v_16(i,1)*(2^7) + bi2de([0 v_16(i,2:8)],'left-msb') + bi2de(v_16(i,9:16),'left-msb')/(2^8);
end

v_tmp_1 = zeros(numel(v),16); v_1 = zeros(numel(v),16); v_2 = zeros(numel(v),16); v_lsli = zeros(numel(v),16);

v_tmp_1(1,:) = zeros(1,16);
v_1(1,:) = bin_add(sgn_xtnd(v_16(1,:)),zeros(1,16));
v_2(1,:) = bin_add(v_1(1,:),zeros(1,16));
v_lsli(1,:) = bin_add(v_2(1,:),v_tmp_1(1,:));

tmp1 = bin_add(v_tmp16(1,:),two_cmpl(zeros(1,16)));
tmp2 = bin_add(rgt_xtnd(tmp1),tmp1);
tmp3 = bin_add(rgt_xtnd(rgt_xtnd(rgt_xtnd(tmp1))),tmp2);
tmp4 = bin_add(rgt_xtnd(rgt_xtnd(rgt_xtnd(rgt_xtnd(tmp1)))),tmp3);
v_tmp_1(2,:) = bin_add(rgt_xtnd(rgt_xtnd(rgt_xtnd(rgt_xtnd(rgt_xtnd(tmp1))))),tmp4);

v_1(2,:) = bin_add(sgn_xtnd(v_16(2,:)),two_cmpl(sgn_xtnd(v_16(1,:))));
v_2(2,:) = bin_add(v_1(2,:),two_cmpl(v_1(1,:)));
v_lsli(2,:) = bin_add(v_2(2,:),v_tmp_1(2,:));

for i=3:numel(v)
tmp1 = bin_add(v_tmp16(i-1,:),two_cmpl(v_tmp16(i-2,:)));
tmp2 = bin_add(rgt_xtnd(tmp1),tmp1);
tmp3 = bin_add(rgt_xtnd(rgt_xtnd(rgt_xtnd(tmp1))),tmp2);
tmp4 = bin_add(rgt_xtnd(rgt_xtnd(rgt_xtnd(rgt_xtnd(tmp1)))),tmp3);
v_tmp_1(i,:) = bin_add(rgt_xtnd(rgt_xtnd(rgt_xtnd(rgt_xtnd(rgt_xtnd(tmp1))))),tmp4);
    
v_1(i,:) = bin_add(sgn_xtnd(v_16(i,:)),two_cmpl(sgn_xtnd(v_16(i-1,:))));
v_2(i,:) = bin_add(v_1(i,:),two_cmpl(v_1(i-1,:)));
v_lsli(i,:) = bin_add(v_2(i,:),v_tmp_1(i,:));
end

v_tmp_1_a = zeros(1,numel(v));
for i=1:numel(v)
v_tmp_1_a(i) = -v_tmp_1(i,1)*(2^7) + bi2de([0 v_tmp_1(i,2:8)],'left-msb') + bi2de(v_tmp_1(i,9:16),'left-msb')/(2^8);
end
v_1_a = zeros(1,numel(v));
for i=1:numel(v)
v_1_a(i) = -v_1(i,1)*(2^7) + bi2de([0 v_1(i,2:8)],'left-msb') + bi2de(v_1(i,9:16),'left-msb')/(2^8);
end
v_2_a = zeros(1,numel(v));
for i=1:numel(v)
v_2_a(i) = -v_2(i,1)*(2^7) + bi2de([0 v_2(i,2:8)],'left-msb') + bi2de(v_2(i,9:16),'left-msb')/(2^8);
end
v_lsli_a = zeros(1,numel(v));
for i=1:numel(v)
v_lsli_a(i) = -v_lsli(i,1)*(2^7) + bi2de([0 v_lsli(i,2:8)],'left-msb') + bi2de(v_lsli(i,9:16),'left-msb')/(2^8);
end

spec = fft(v_tmp_a.*ds_hann(N))/(N/4);
figure;
subplot(2,1,1);
plot(log10(f),dbv(spec(1:N/2+1))); hold on;
snr = calculateSNR(spec(3:fB+5),fB-2);
NBW = 1.5/N;
Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / 3;
plot(log10(f),dbp(Sqq*NBW),'m','Linewidth',2); grid on;
text(0.001, -90, sprintf('NBW = %4.1E x f_s ',NBW),'Hor','right');
legend(sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),'1-bit Expected');

subplot(2,1,2);
spec_v = fft(v_a.*ds_hann(N))/(N/4);
plot(log10(f),dbv(spec_v(1:N/2+1))); hold on;
snr_v = calculateSNR(spec_v(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
plot(log10(f),dbv(spec(1:N/2+1))); hold on;
plot(log10(f),dbv(spec_lsli(1:N/2+1))); grid on;
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);

legend(sprintf('2-bit Quantizer SNR = %4.1fdB @ OSR = %d',snr_v,OSR),sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),sprintf('Leslie Singh SNR = %4.1fdB @ OSR = %d',snr_lsli,OSR));
