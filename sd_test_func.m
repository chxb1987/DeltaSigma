function [snrv,v,y,yt,vt] = sd_test_func (N,fs,f,u,w,order,OSR,ntf,nlev) 
ut = u+w;
[vw,xnw,xmaxw,yw] = simulateDSM(w,ntf,nlev);
[v,xn,xmax,y] = simulateDSM(u,ntf,nlev);
[vt,xnt,xmaxt,yt] = simulateDSM(ut,ntf,nlev);
fx = -0.5*fs:fs/N:0.5*fs; fint = (0.5*N+3):((0.5*N+1)+N*1.5*f/fs); 
hwv = fftshift(abs(fft(v'.*hanning(N+1))))./(N+1); hwvw = fftshift(abs(fft(vw'.*hanning(N+1))))./(N+1);
hwut = fftshift(abs(fft(ut'.*hanning(N+1))))./(N+1); hwvt = fftshift(abs(fft(vt'.*hanning(N+1))))./(N+1);
hwu = fftshift(abs(fft(u'.*hanning(N+1))))./(N+1); hwuw = fftshift(abs(fft(w'.*hanning(N+1))))./(N+1);
snru = 10*log10(sum(hwu(fint).^2)/sum(hwuw(fint).^2)); snrv = 10*log10(sum(hwv(fint).^2)/sum(hwvw(fint).^2));
plot(log10(fx),20.*log10(hwut)); hold on; plot(log10(fx),20.*log10(hwvt)); grid on;
legend(sprintf('snru = %f',snru),sprintf('snrv = %f',snrv));
end
