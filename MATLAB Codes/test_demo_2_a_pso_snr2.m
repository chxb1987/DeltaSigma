tic;
n = 50; c1 = 0.5; c2 = 0.5; w = linspace(1,0,200); itrmax = 200; gmin = 1e-100; 
N = 32768; OSR = N/6; nlev = 2;
z = zpk('z',1);

H = ((z-1)^2)/(z^2-1.567*z+0.7835);

% gamma = 1e3;
% theta = 6;
% delta = 1;
% alpha = gamma/(1+delta*(1+gamma));
% beta = delta*(1+gamma)/(1+delta*(1+gamma));
% H = ((z-beta)^2)/(z^2-(2*beta-0.433*alpha)*z+(beta^2+0.2165*alpha^2-0.433*alpha*beta));

fB = ceil(N/(2*OSR)); 
u_theo = 0.586*sin(2*pi*fB/N*(0:N-1));	% half-scale sine-wave input
% [vtheo,~,~,y_theo] = simulateDSM(u_theo,H,nlev); 

% sdin = csvread('JB_Q.csv',1,0);
% sdin(:,2) = 2.*(double(sdin(:,2) > 0.9))-1;
% vtheo = sdin(100:(N+99),2)';

% sdin = csvread('Q16384_4M.csv',1,0);
% sdin_m(:,1) = sdin(2:4:size(sdin,1),1);
% sdin_m(:,2) = sdin(2:4:size(sdin,1),2);
% sdin_m(:,2) = 2.*(double(sdin_m(:,2) > 0.9))-1;
% vtheo = sdin_m(100:(N+99),2)';

sdin = csvread('Q_3_4M.csv',1,0);
sdin_m(:,1) = sdin(2:4:size(sdin,1),1);
sdin_m(:,2) = sdin(2:4:size(sdin,1),2);
sdin_m(:,2) = 2.*(double(sdin_m(:,2) > 0.9))-1;
vtheo = sdin_m(100:(N+99),2)';

vtheo_a = vtheo./2;

% yin1 = csvread('JB_y.csv',1,0);
% yin1(:,2) = (yin1(:,2) - 0.8985)./0.1;
% y_theo = yin1(100:(N+99),2)';

% yin = csvread('von16384_4M.csv',1,0);
% yin_m(:,1) = yin(2:4:size(yin,1),1);
% yin_m(:,2) = yin(2:4:size(yin,1),2);
% yin_m(:,2) = (yin_m(:,2) - 0.9140)./0.1;
% y_theo = yin_m(100:(N+99),2)';

yin = csvread('von_3_4M.csv',1,0);
yin_m(:,1) = yin(2:4:size(yin,1),1);
yin_m(:,2) = yin(2:4:size(yin,1),2);
yin_m(:,2) = (yin_m(:,2) - 0.9140)./0.1;
y_theo = yin_m(100:(N+99),2)';

for i=1:n
    x(1,:,i) = [0.4 + 0.2*rand, 0.3 + 0.2*rand, 0.9 + 0.1*rand];
    u(1,:,i) = 0.1.*[0.4 + 0.2*rand, 0.3 + 0.2*rand, 0.9 + 0.1*rand];
end

P(1,:,:) = x(1,:,:); optim = zeros(1,3);
g = ones(itrmax,n); gG = ones(1,itrmax); gpl = zeros(1,itrmax);
for i=1:n
    g(1,i) = pso_gfunc2(x(1,:,i),vtheo_a,y_theo,0);    
end
[~,ind] = min(g(1,:));
G(1,:) = P(1,:,ind);
gpl(1) = min(g(1,:));

for j=2:itrmax
    for i=1:n
        k = 2./abs(2-(c1+c2)-sqrt((c1+c2)^2-4*(c1+c2)));
%         u(j,:,i) = k.*(u(j-1,:,i) + (c1*rand).*(P(j-1,:,i) - x(j-1,:,i)) + (c2*rand).*(G(j-1,:) - x(j-1,:,i)));
        u(j,:,i) = w(j).*u(j-1,:,i) + (c1*rand).*(P(j-1,:,i) - x(j-1,:,i)) + (c2*rand).*(G(j-1,:) - x(j-1,:,i));
        x(j,:,i) = x(j-1,:,i) + u(j,:,i);
        
        x(j,:,i) = x(j,:,i) + [(1e-3).*rand(1,2),0];
        u(j,:,i) = u(j,:,i) + [(1e-4).*rand(1,2),0];
        
        g(j,i) = pso_gfunc2(x(j,:,i),vtheo_a,y_theo,0);
        if (g(j,i) < g(j-1,i))
            P(j,:,i) = x(j,:,i);
        else
            P(j,:,i) = P(j-1,:,i);
            g(j,i) = g(j-1,i);
        end
    end
    [~,ind] = min(g(j,:));
    G(j,:) = P(j,:,ind);
    
    gG(j) = pso_gfunc2(G(j,:),vtheo_a,y_theo,0);
    if ((j>3)&&(j~=itrmax))
        w(j+1) = exp(-mean(abs(g(j-1,:)-gG(j-1)))./mean(abs(g(j-2,:)-gG(j-2))));
%         w(j+1) = 1./(exp(mean(abs(g(j,:)-gG(j)))-mean(abs(g(j-1,:)-gG(j-1))))+1);
        if (isnan(w(j+1)))
            w(j+1) = 0;
        end
    end
    if (min(g(j,:)) < gmin)
        optim = G(j,:);
        break;
    end
    gpl(j) = min(g(j,:));
end

if (sum(optim == zeros(1,3)) == 3)
    optim = G(itrmax,:);
end
mmg = min(min(g));

% figure;
% plot(log10(gpl),'k'); grid on; hold on;
% title('Cost Function vs Number of Iterations')
% xlabel('Number of Iterations')
% ylabel('log10(g)')

% figure;

% pso_gfunc_lsli([1,1.567,0.7835],vtheo_a,y_theo,1); hold on;
% pso_gfunc_lsli([beta,(2*beta-0.433*alpha),(beta^2+0.2165*alpha^2-0.433*alpha*beta)],vtheo_a,y_theo,1); hold on;

pso_gfunc2(optim,vtheo_a,y_theo,1);
title('SNR Improvement vs Number of Quantizer Bits for Leslie-Singh Topology')
xlabel('B (Quantizer Bits)')
ylabel('SNR Improvement (dB)')
% legend('Known Coefficients','Coefficients Obtained from PSO')
toc;
                
