N = 32768; fs = 1e6; f = (33*fs)/N; Au = [-100:5:-10 -9:1:0]; snr = zeros(1,numel(Au));
u = 3*sin(2*pi*f.*(0:1/fs:N/fs)); w = (10e-4).*normrnd(0,1,1,N+1); 
order = 2; OSR = N/66; nlev = 2;
ntf = synthesizeNTF(order,OSR); plotPZ(ntf);
for i=1:numel(Au)
    [snr(i),v] = sd_test_func(N,fs,f,10^(Au(i)/20)*u,w,order,OSR,ntf,nlev);    
end
[sqnr, amp] = simulateSNR(ntf,OSR,[],[],nlev);
figure;
plot(Au,snr); hold on; plot(amp,sqnr); grid on;
form = 'CIFB'; [a,g,b,c] = realizeNTF(ntf,form); b(2:end) = 0;
ABCD = stuffABCD(a,g,b,c,form);
ABCDs = scaleABCD(ABCD,nlev,0,0.5,[],0.9);
[a,g,b,c] = mapABCD(ABCDs,form);


