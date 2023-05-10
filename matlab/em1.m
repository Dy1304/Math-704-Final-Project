%EM Euler-Maruyama method on linear SDE % % SDE is dX = mu*X dt + sigma*X dW, X(0) = Xzero, % where mu = 2, sigma = 1 and Xzero = 1.

% % Discretized Brownian path over [0,1] has dt = 2^(-8). % Euler-Maruyama uses timestep R*dt.
clear all; clc;

AAPL1 = readmatrix('AAPL.csv');
AAPL = AAPL1(:,2);
Xzero = AAPL(1);
n = length(AAPL); Dt = 1/n; T = 1;

randn('state',100)
mu = 0.1244; sigma = 0.1038; % problem parameters
N = 2^12; dt = 1/N; % Brownian increments 
dW = sqrt(dt)*randn(1,N); W = cumsum(dW); % discretized Brownian path

plot([Dt:Dt:T], AAPL,'g-'), hold on % plot the true stock price

Xtrue = Xzero*exp((mu-0.5*sigma^2)*([dt:dt:T])+sigma*W); 
plot([0:dt:T],[Xzero,Xtrue],'m-'), hold on % plot the true solution 

emerr1 = ones(1,6);
for i = 1:6
 R = 2^(i);
 Dt = R*dt; L = N/R; % L EM steps of size Dt = R*dt

Xem = zeros(1,L); Xtemp = Xzero; % preallocate for efficiency
for j = 1:L 
Winc = sum(dW(R*(j-1)+1:R*j));
Xtemp = Xtemp + Dt*mu*Xtemp + sigma*Xtemp*Winc;
Xem(j) = Xtemp; 
end

emerr = abs(Xem(end)-Xtrue(end));
emerr1(i) = emerr;
end

plot([0:Dt:T],[Xzero,Xem],'r--*'), hold off 
xlabel('t','FontSize',12) 
ylabel('X','FontSize',16,'Rotation',0,'HorizontalAlignment','right')
legend('real stock price','true solution','EM method solution')
emerr1