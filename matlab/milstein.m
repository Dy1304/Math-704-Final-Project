clear all; clc;
AAPL1 = readmatrix('AAPL.csv');
AAPL = AAPL1(:,2); % read data
Xzero = AAPL(1);

rand('state',100) 
%r = 2; K = 1; beta = 0.25; Xzero = 0.5; % problem parameters
mu = 0.1244; sigma = 0.1038; % problem parameters
T = 1; N = 2^(11); dt = T/N; 
M = 500; R = [1; 16; 32; 64; 128]; % number of paths sampled

dW = sqrt(dt)*randn(M,N); % Brownian increments
Xmil = zeros(M,5); 
for p = 1:5
 Dt = R(p)*dt;  % L timesteps of size Dt = R dt
 L = N/R(p); 
 Xtemp = Xzero*ones(M,1); 
    for j = 1:L
        Winc = sum(dW(:,R(p)*(j-1)+1:R(p)*j),2);
        Xtemp = Xtemp + mu*Xtemp*Dt + sigma*Xtemp.*Winc+ 0.5*sigma^2*Xtemp.*(Winc.^2 - Dt); 
    end
    Xmil(:,p) = Xtemp; % store Milstein solution at t =1
end

Xref = Xmil(:,1); % Reference solution
Xerr = abs(Xmil(:,2:5) - repmat(Xref,1,4)); % Error in each path
mean(Xerr); % Mean pathwise erorrs
Dtvals = dt*R(2:5); % Milstein timesteps used

loglog(Dtvals,mean(Xerr),'b*-'), hold on 
loglog(Dtvals,Dtvals,'r--'), hold off % reference slope of 1
axis([1e-3 1e-1 1e-4 1]) 
xlabel('\Delta t')
ylabel('Sample average of | X(T) - X_L |') 
title('milstein2','FontSize',10)

%%%% Least squares fit of error = C * Dt^q %%%% 
A = [ones(4,1), log(Dtvals)]; 
rhs = log(mean(Xerr)'); 
sol = A\rhs; q = sol(2) 
resid = norm(A*sol - rhs)