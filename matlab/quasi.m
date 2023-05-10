clear all;

syms t d
S = sym(10);        % current stock price (spot price)
K = sym(12);         % exercise price (strike price)
sigma = sym(0.50);   % volatility of stock
T = sym(0.5);       % expiry time in years
r = sym(0.05);       % annualized risk-free interest rate

PV_K = K*exp(-r*T);
d1 = (log(S/K) + (r + sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
Nd1 = int(exp(-((t)^2)/2),-Inf,d1)*1/sqrt(2*pi);
Nd2 = int(exp(-((t)^2)/2),-Inf,d2)*1/sqrt(2*pi);
C = Nd1*S - Nd2*PV_K;

Csym = Nd1*S - Nd2*PV_K;
Cvpa = vpa(Csym)


[Call, Put] = blsprice(double(S), double(K), double(r), double(T), double(sigma))

rand('seed',10);
U_1 = rand(1,100000);
U_2 = rand(1,100000);
R0 = sqrt(-2*log(U_1));
theta = 2*pi*U_2;
x = R0.*cos(theta);
y = R0.*sin(theta); % generate quasi-random numbers

randn('state', 100)
mu = 0.05; sigma = 0.5; Xzero = 10; r = 0.05; K = 12; % problem parameters
T = 0.5; N = 2^11; dt = T/N;
M = 10000; k = 5; % number of sample paths and test number of intervals

Xerr1 = zeros(M,k); 
for s = 1:M
    dW = sqrt(dt)* randn(1,N); % pseudo Brownian increments
    W = cumsum(dW); % pseudo discrete Brownian path
    dQ = sqrt(dt)* y(1:N); % quasi Brownian increments
    Q = cumsum(dQ); % quasi discrete Brownian path
    Xtrue = 0.81722509732109537726399447901733;
    for p = 1:k 
        R = 2^(p-1); Dt = R*dt; L = N/R; % L Euler steps of size Dt = R*dt 
        Xtemp = Xzero; Xtemp1 = Xzero; 
        Xval1 = zeros(size(Xtemp)); Xval2 = zeros(size(Xtemp));
        Winc = sqrt(Dt)* randn(1,L);
        WinQ = sqrt(Dt)* y(1:L);
        for j = 1:L 
            %Winc = sum(dW(R*(j-1)+1:R*j)); 
            %WinQ = sum(dQ(R*(j-1)+1:R*j));
            Xtemp = Xtemp + Dt*mu*Xtemp + sigma*Xtemp*Winc(j); % EM method
            Xtemp1 = Xtemp1 + Dt*mu*Xtemp1 + sigma*Xtemp1*WinQ(j); % quasi-random number
            Xval1(j) = Xtemp; Xval2(j) = Xtemp1;
        end
        X1 = exp(-r*T)*mean(max((Xval1-K),0)); % mean(max((Xval1-K),0)) for European call option
        X2 = exp(-r*T)*mean(max((Xval2-K),0));
        Xerr1(s,p) = abs(X1 - Xtrue); 
        Xerr2(s,p) = abs(X2 - Xtrue); 
    end
end

mean(Xerr1)
mean(Xerr2)

Dtvals = dt*(2.^([0:4]));
loglog(Dtvals,mean(Xerr1),'b*-'), hold on
loglog(Dtvals,mean(Xerr2),'ro-'), hold off % reference slope of 1/2
axis([1e-4 1e-2 1e-2 1.2])
xlabel('\Delta t'), ylabel('Error')
% title('strongcon_em.m','FontSize',10)
title('error produced by pseudo and quasi-random numbers','FontSize',10)

