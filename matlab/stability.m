clear all;

randn('state',100) 
T = 200; M = 50000; Xzero = 1; % set T, number of paths and initial value
ltype = {'b-','r--','m-.'};

%%%%%%%%%%%% Mean Square %%%%%%%%%%%%%
subplot(211)  
mu = 0.1244; sigma = 0.1038; % problem parameters 
for k = 1:3
    Dt = 2^(1-k); 
    N = T/Dt; 
    Xms = zeros(1,N); 
    Xtemp = Xzero*ones(M,1); 
    for j = 1:N
        Winc = sqrt(Dt)*randn(M,1);
        Xtemp = Xtemp + Dt*mu*Xtemp + sigma*Xtemp.*Winc;
        Xms(j) = mean(Xtemp.^2); % mean-square estimate 
    end
    semilogy([0:Dt:T],[Xzero,Xms],ltype{k},'Linewidth',2), hold on
end
legend('\Delta t = 1','\Delta t = 1/2','\Delta t = 1/4') 
title('Mean-Square: \mu = -1, \sigma = \surd 0.5','FontSize',16) 
ylabel('E[X^2]','FontSize',12), axis([0,T,1e-20,1e+20]), hold off

subplot(212) %%%%% Asymptotic: a single path %%%%%%% 
T = 500; mu = 0.1244; sigma = 0.1038; % problem parameters 
for k = 1:3
    Dt = 2^(1-k); 
    N = T/Dt; 
    Xemabs = zeros(1,N); Xtemp = Xzero; 
    for j = 1:N
        Winc = sqrt(Dt)*randn;
        Xtemp = Xtemp + Dt*mu*Xtemp + sigma*Xtemp*Winc;
        Xemabs(j) = abs(Xtemp); 
    end
    semilogy([0:Dt:T],[Xzero,Xemabs],ltype{k},'Linewidth',2), hold on
end
legend('\Delta t = 1','\Delta t = 1/2','\Delta t = 1/4') 
title('Single Path: \mu = 0.02, \sigma = \surd 2','FontSize',16) 
ylabel('|X|','FontSize',12), axis([0,T,1e-50,1e+100]), hold off