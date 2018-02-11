function g = pso_gfunc (optim,vtheo_a,v_a)
N = 32768; OSR = N/66; fB = ceil(N/(2*OSR)); 

v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];

v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;

spec = fft(vtheo_a.*ds_hann(N))/(N/4);
snr = calculateSNR(spec(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);

g = 10^((-snr_lsli+snr));
end