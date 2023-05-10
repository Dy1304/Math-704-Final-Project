clear all; clc;
AAPL1 = readmatrix('AAPL.csv');
AAPL = AAPL1(:,2); % read data
Xzero = AAPL(1);

randn('state',100); 
mu = 2; sigma = 0.1; Xzero = 1;
%mu = 0.1244; sigma = 0.1038;  problem parameters
T = 1; k = 5;
M = 50000; % number of paths sampled

Xem = zeros(k,1); % preallocate arrays
for p = 1:k % take various Euler timesteps
    Dt = 2^(p-11); L = T/Dt; % L Euler steps of size Dt
    Xtemp = Xzero*ones(M,1); 
        for j = 1:L
            Winc = sqrt(Dt)*randn(M,1);
            Winc = sqrt(Dt)*sign(randn(M,1)); %% use for weak E-M %%
            Xtemp = Xtemp + Dt*mu*Xtemp + sigma*Xtemp.*Winc; 
        end
        Xem(p) = mean(Xtemp);
end
Xerr = abs(Xem - Xzero*exp(mu));

Dtvals = 2.^([1:k]-11);  
loglog(Dtvals,Xerr,'b*-'), hold on 
loglog(Dtvals,Dtvals,'r--'), hold off % reference slope of 1 
axis([1e-3 1e-1 1e-4 1]) 
xlabel('\Delta t'), ylabel('| E(X(T)) - Sample average of X_L |') 
title('emweak.m','FontSize',10)

%%%% Least squares fit of error = C * dt^q %%%% 
A = [ones(p,1), log(Dtvals)']; rhs = log(Xerr); 
sol = A\rhs; 
q = sol(2) 
resid = norm(A*sol - rhs)
