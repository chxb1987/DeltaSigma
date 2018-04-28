function g = pso_gfunc2 (prm,vtheo_a,y_theo,flag)
N = 32768; OSR = N/6; fB = ceil(N/(2*OSR)); 
optim = [prm(3),2*prm(3)-prm(2),prm(1)*prm(2)-prm(3)*prm(2)+prm(3)^2];

v = ds_quantize(((2^4)-1).*y_theo,2^4); 
v_a = v./((2^4));    
% v_a = y_theo;
    
v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];

v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

v_lsli_a = v_tmp_1_a + v_2_a;

spec = fft(vtheo_a.*ds_hann(N))/(N/4);
snr = calculateSNR(spec(3:fB+5),fB-2);

spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);

g = (10^((-snr_lsli+snr)))./(1+10^((-snr_lsli+snr)));
% g = 10^((-snr_lsli+snr));
if (flag == 1)
    interm = zeros(1,6);
    for i=1:6
    v = ds_quantize(((2^(i+1))-1).*y_theo,2^(i+1)); 
    v_a = v./(2^(i+1));
    v_tmp_1_a = (2.*optim(1)-optim(2)).*[0 vtheo_a(1:numel(v_a)-1)] - (optim(1)^2-optim(3)).*[0 0 vtheo_a(1:numel(v_a)-2)];
    
    v_1_a = 0.5.*(v_a - optim(1).*[0 v_a(1:numel(v_a)-1)]);
    v_2_a = (v_1_a - optim(1).*[0 v_1_a(1:numel(v_a)-1)]);

    v_lsli_a = v_tmp_1_a + v_2_a;

    spec = fft(vtheo_a.*ds_hann(N))/(N/4);
    snr = calculateSNR(spec(3:fB+5),fB-2);

    spec_lsli = fft(v_lsli_a.*ds_hann(N))/(N/4);
    snr_lsli = calculateSNR(spec_lsli(3:fB+5),fB-2);
    interm(i) = snr_lsli-snr;
    end
plot(2:7,interm,'--k'); grid on; 
end
end