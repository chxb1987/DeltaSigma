function [optim,w,mmg] = test_demo_2_a_pso 
tic;
n = 20; c1 = 1.5; c2 = 1.5; w = 0.5+0.5.*rand(1,200); itrmax = 200; gmin = 0.01; 
N = 32768; OSR = N/66; nlev = 2;
z = zpk('z',1);
fB = ceil(N/(2*OSR)); ftest = floor(2/3*fB); f = linspace(0,0.5,N/2+1);
y = 0.7071*sin(2*pi*fB/N*(0:N-1));	% half-scale sine-wave input
% a = [100 -100 -100 2.5 0.75 0.01 3 -0.1 50 0.5 0.05 -100 25 -0.1 100 80 0.02 -200 1 0.02 -200 2 -0.2 -100];
% ptheo = [55.17,-72.14,-49.42,1.2,0.36,0.003,1.5,-0.03,25,0.22,0.018,-52,15.75,-0.03,45,40,0.075,-90,0.57,0.065,-90,1,-0.1,-30];

sdin = csvread('Q_op3.csv',1,0);
sdin(:,2) = 2.*(double(sdin(:,2) > 0.9))-1;
vtheo = sdin(10:(N+9),2)';
% sdin = SBBoser;
% sdin = 2.*(double(sdin > 0.9))-1;
% vtheo = sdin(10:(N+9))';
% H1 = ((z-0.96)^2)/(z^2-1.433*z+0.7585);
% vtheo = simulateDSM(y,H1,nlev);
spec_theo = fft(vtheo.*ds_hann(N))/(N/4);

yin1 = csvread('von_op3.csv',1,0);
yin2 = csvread('vip1_op3.csv',1,0);
yin1(:,2) = (yin1(:,2) - 0.9)./0.1;
yin2(:,2) = (yin2(:,2) - 0.9)./0.125;
% spec_theo = fft(yin1(10:(N+9),2)'.*ds_hann(N))/(N/4);

% x(1,:,:) = [100*rand(1,1,n) -100*rand(1,1,n) -100*rand(1,1,n) 2.5*rand(1,1,n) 0.75*rand(1,1,n) 0.01*rand(1,1,n) 3*rand(1,1,n) -0.1*rand(1,1,n) 50*rand(1,1,n) 0.5*rand(1,1,n) 0.05*rand(1,1,n) -100*rand(1,1,n) 25*rand(1,1,n) -0.1*rand(1,1,n) 100*rand(1,1,n) 80*rand(1,1,n) 0.02*rand(1,1,n) -200*rand(1,1,n) rand(1,1,n) 0.02*rand(1,1,n) -200*rand(1,1,n) 2*rand(1,1,n) -0.2*rand(1,1,n) -100*rand(1,1,n)];
% u(1,:,:) = 0.1.*[100*rand(1,1,n) -100*rand(1,1,n) -100*rand(1,1,n) 2.5*rand(1,1,n) 0.75*rand(1,1,n) 0.01*rand(1,1,n) 3*rand(1,1,n) -0.1*rand(1,1,n) 50*rand(1,1,n) 0.5*rand(1,1,n) 0.05*rand(1,1,n) -100*rand(1,1,n) 25*rand(1,1,n) -0.1*rand(1,1,n) 100*rand(1,1,n) 80*rand(1,1,n) 0.02*rand(1,1,n) -200*rand(1,1,n) rand(1,1,n) 0.02*rand(1,1,n) -200*rand(1,1,n) 2*rand(1,1,n) -0.2*rand(1,1,n) -100*rand(1,1,n)];

for i=1:n
%     x(1,:,i) = [0.95 + 0.05*rand, -sqrt(2) + 2*sqrt(2)*rand, 0.5*rand];
%     u(1,:,i) = 0.1.*[0.95 + 0.05*rand, -sqrt(2) + 2*sqrt(2)*rand, 0.5*rand];
    x(1,:,i) = [0.95 + 0.05*rand, 1.5 + 0.1*rand, 0.75 + 0.1*rand];
    u(1,:,i) = 0.1.*[0.95 + 0.05*rand, 1.5 + 0.1*rand, 0.75 + 0.1*rand];
end

P(1,:,:) = x(1,:,:); optim = zeros(1,3);
g = ones(itrmax,n); gG = ones(1,itrmax);
for i=1:n
    H1 = ((z-x(1,1,i))^2)/(z^2-x(1,2,i)*z+x(1,3,i));
%     H1 = ((z-x(1,1,i))^2)/(z^2-1.5*z+0.75);        
    [tmp,~,~,tmp_y] = simulateDSM(y,H1,nlev); 
    spec = fft(tmp.*ds_hann(N))/(N/4);    
    g(1,i) = rms(dbv(spec(1:N/2+1))-dbv(spec_theo(1:N/2+1)))/rms(dbv(spec_theo(1:N/2+1)));    
%     tmp = hh_rk4_script_new(-60,I,x(1,:,i));    
%     g(1,i) = sqrt(mean([(tmp(1:80) - vtheo(1:80)).*(tmp(1:80) - vtheo(1:80));2.*(tmp(81:139) - vtheo(81:139)).*(tmp(81:139) - vtheo(81:139));(tmp(140:1001) - vtheo(140:1001)).*(tmp(140:1001) - vtheo(140:1001))]))/sqrt(mean([vtheo(1:80).*vtheo(1:80);2.*vtheo(81:139).*vtheo(81:139);vtheo(140:1001).*vtheo(140:1001)]));
%     g(1,i) = rms(tmp-vtheo)/rms(vtheo);
    if (isnan(g(1,i)))
       g(1,i) = 1;
    end
end
[~,ind] = min(g(1,:));
G(1,:) = P(1,:,ind);

for j=2:itrmax
    for i=1:n
        u(j,:,i) = w(j).*u(j-1,:,i) + (c1*rand).*(P(j-1,:,i) - x(j-1,:,i)) + (c2*rand).*(G(j-1,:) - x(j-1,:,i));
        x(j,:,i) = x(j-1,:,i) + u(j,:,i);
        
%         tmpx = x(j,:,i).*a - x(j,:,i).*x(j,:,i);
%         if (sum(logical(tmpx < 0)) > 0)
%             u(j,logical(tmpx < 0),i) = (0.1.*a(logical(tmpx < 0))).*rand(1,sum(logical(tmpx < 0)));
%             x(j,logical(tmpx < 0),i) = a(logical(tmpx < 0)).*rand(1,sum(logical(tmpx < 0)));            
%         end
%         if ((abs(-x(1,1,i) + sqrt(x(1,1,i)^2 - 4*x(1,2,i))) >= 2)||(abs(-x(1,1,i) - sqrt(x(1,1,i)^2 - 4*x(1,2,i))) >= 2))
%             x(j,:,i) = x(j-1,:,i);
%         end
        Hj = ((z-x(j,1,i))^2)/(z^2-x(j,2,i)*z+x(j,3,i));
%         Hj = ((z-x(j,1,i))^2)/(z^2-1.5*z+0.75);
        [tmp,~,~,tmp_y] = simulateDSM(y,Hj,nlev);
        spec = fft(tmp.*ds_hann(N))/(N/4);
        g(j,i) = rms(dbv(spec(1:N/2+1))-dbv(spec_theo(1:N/2+1)))/rms(dbv(spec_theo(1:N/2+1)));
%         tmp = hh_rk4_script_new(-60,I,x(j,:,i));
%         g(j,i) = sqrt(mean([(tmp(1:80) - vtheo(1:80)).*(tmp(1:80) - vtheo(1:80));2.*(tmp(81:139) - vtheo(81:139)).*(tmp(81:139) - vtheo(81:139));(tmp(140:1001) - vtheo(140:1001)).*(tmp(140:1001) - vtheo(140:1001))]))/rms(vtheo);
%         g(j,i) = sqrt(mean([(tmp(1:80) - vtheo(1:80)).*(tmp(1:80) - vtheo(1:80));2.*(tmp(81:139) - vtheo(81:139)).*(tmp(81:139) - vtheo(81:139));(tmp(140:1001) - vtheo(140:1001)).*(tmp(140:1001) - vtheo(140:1001))]))/sqrt(mean([vtheo(1:80).*vtheo(1:80);2.*vtheo(81:139).*vtheo(81:139);vtheo(140:1001).*vtheo(140:1001)]));
%         g(j,i) = rms(tmp-vtheo)/rms(vtheo);
        if (isnan(g(j,i)))
            g(j,i) = 1;
        end
        if (g(j,i) < g(j-1,i))
            P(j,:,i) = x(j,:,i);
        else
            P(j,:,i) = P(j-1,:,i);
            g(j,i) = g(j-1,i);
        end
    end
    [~,ind] = min(g(j,:));
    G(j,:) = P(j,:,ind);
    HG = ((z-G(j,1))^2)/(z^2-G(j,2)*z+G(j,3));
%     HG = ((z-G(j,1))^2)/(z^2-1.5*z+0.75);
    [tmpG,~,~,tmpG_y] = simulateDSM(y,HG,nlev);
    specG = fft(tmpG.*ds_hann(N))/(N/4);
    gG(j) = rms(dbv(specG(1:N/2+1))-dbv(spec_theo(1:N/2+1)))/rms(dbv(spec_theo(1:N/2+1)));
%     tmpG = hh_rk4_script_new(-60,I,G(j,:));
%     gG(j) = sqrt(mean([(tmpG(1:80) - vtheo(1:80)).*(tmpG(1:80) - vtheo(1:80));2.*(tmpG(81:139) - vtheo(81:139)).*(tmpG(81:139) - vtheo(81:139));(tmpG(140:1001) - vtheo(140:1001)).*(tmpG(140:1001) - vtheo(140:1001))]))/rms(vtheo);
%     gG(j) = sqrt(mean([(tmpG(1:80) - vtheo(1:80)).*(tmpG(1:80) - vtheo(1:80));2.*(tmpG(81:139) - vtheo(81:139)).*(tmpG(81:139) - vtheo(81:139));(tmpG(140:1001) - vtheo(140:1001)).*(tmpG(140:1001) - vtheo(140:1001))]))/sqrt(mean([vtheo(1:80).*vtheo(1:80);2.*vtheo(81:139).*vtheo(81:139);vtheo(140:1001).*vtheo(140:1001)]));
%     gG(j) = rms(tmpG-vtheo)/rms(vtheo);
    if (j>3)
        w(j+1) = exp(-mean(abs(g(j-1,:)-gG(j-1)))./mean(abs(g(j-2,:)-gG(j-2))));
        if (isnan(w(j+1)))
            w(j+1) = 0;
        end
    end
    if (min(g(j,:)) < gmin)
        optim = G(j,:);
        break;
    end
end

if (sum(optim == zeros(1,3)) == 3)
    optim = G(itrmax,:);
end
mmg = min(min(g));
toc;
end
                
% N = 32768; OSR = N/66; nlev = 2;
% z = zpk('z',1);
% H = (z^2-2*z+1)/(z^2-1.225*z+0.4415); 
% fB = ceil(N/(2*OSR)); ftest = floor(2/3*fB);
% u = 0.5*sin(2*pi*ftest/N*(0:N-1));	% half-scale sine-wave input
% [v_tmp,xn,xnmax,y] = simulateDSM(u,H,nlev); 
% f = linspace(0,0.5,N/2+1);
% spec = fft(v_tmp.*ds_hann(N))/(N/4);
% figure;
% plot(log10(f),dbv(spec(1:N/2+1))); hold on;
% snr = calculateSNR(spec(3:fB+1),ftest-2);
% NBW = 1.5/N;
% Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / 3;
% plot(log10(f),dbp(Sqq*NBW),'m','Linewidth',2); grid on;
% text(0.001, -90, sprintf('NBW = %4.1E x f_s ',NBW),'Hor','right');
% legend(sprintf('1-bit Simulation SNR = %4.1fdB @ OSR = %d',snr,OSR),'1-bit Expected');

% N = 32768; OSR = N/66; nlev = 2;
% z = zpk('z',1);
% fB = ceil(N/(2*OSR)); ftest = floor(2/3*fB); f = linspace(0,0.5,N/2+1);
% y = 0.5*sin(2*pi*fB/N*(0:N-1));
% 
% H1 = ((z-0.96)^2)/(z^2-1.433*z+0.7585);
% vtheo = simulateDSM(y,H1,nlev);
% sp1 = fft(vtheo.*ds_hann(N))/(N/4);
% plot(log10(f),dbv(sp1(1:N/2+1))); hold on;
% snr1 = calculateSNR(sp1(3:fB+5),fB-2);
% 
% H2 = ((z-optim(1))^2)/(z^2-optim(2)*z+optim(3));
% v = simulateDSM(y,H2,nlev);
% sp2 = fft(v.*ds_hann(N))/(N/4);
% plot(log10(f),dbv(sp2(1:N/2+1))); grid on;
% snr2 = calculateSNR(sp2(3:fB+5),fB-2);
% legend(sprintf('%d',snr1),sprintf('%d',snr2));
