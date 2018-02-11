function [optim,w,mmg] = test_demo_2_a_pso_snr 
tic;
n = 50; c1 = 1.5; c2 = 1.5; w = 0.5+0.5.*rand(1,500); itrmax = 500; gmin = 1e-6; 
N = 32768; 

% OSR = N/66; nlev = 2;
% z = zpk('z',1);
% H = ((z-1)^2)/(z^2-1.5*z+0.75);
% fB = ceil(N/(2*OSR)); 
% u_theo = 0.7071*sin(2*pi*fB/N*(0:N-1));	% half-scale sine-wave input
% [vtheo,~,~,y_theo] = simulateDSM(u_theo,H,nlev); 

sdin = csvread('Q_op3.csv',1,0);
sdin(:,2) = 2.*(double(sdin(:,2) > 0.9))-1;
vtheo = sdin(10:(N+9),2)';
vtheo_a = vtheo./2;

yin1 = csvread('von_op3.csv',1,0);
yin1(:,2) = (yin1(:,2) - 0.9)./0.1;
y_theo = yin1(10:(N+9),2)';
v = ds_quantize(((2^5)-1).*y_theo,2^5); 
v_a = v./((2^5));

for i=1:n
    x(1,:,i) = [1 + 0.05*randn, 1.5 + 0.1*randn, 0.75 + 0.1*randn];
    u(1,:,i) = 0.1.*[0.95 + 0.05*randn, 1.5 + 0.1*randn, 0.75 + 0.1*randn];
end

P(1,:,:) = x(1,:,:); optim = zeros(1,3);
g = ones(itrmax,n); gG = ones(1,itrmax);
for i=1:n
    g(1,i) = pso_gfunc(x(1,:,i),vtheo_a,v_a);    
end
[~,ind] = min(g(1,:));
G(1,:) = P(1,:,ind);

for j=2:itrmax
    for i=1:n
        u(j,:,i) = w(j).*u(j-1,:,i) + (c1*rand).*(P(j-1,:,i) - x(j-1,:,i)) + (c2*rand).*(G(j-1,:) - x(j-1,:,i));
        x(j,:,i) = x(j-1,:,i) + u(j,:,i);
        
        g(j,i) = pso_gfunc(x(j,:,i),vtheo_a,v_a);
        if (g(j,i) < g(j-1,i))
            P(j,:,i) = x(j,:,i);
        else
            P(j,:,i) = P(j-1,:,i);
            g(j,i) = g(j-1,i);
        end
    end
    [~,ind] = min(g(j,:));
    G(j,:) = P(j,:,ind);
    
    gG(j) = pso_gfunc(G(j,:),vtheo_a,v_a);
%     if (j>3)
%         w(j+1) = exp(-mean(abs(g(j-1,:)-gG(j-1)))./mean(abs(g(j-2,:)-gG(j-2))));
%         if (isnan(w(j+1)))
%             w(j+1) = 0;
%         end
%     end
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
                
