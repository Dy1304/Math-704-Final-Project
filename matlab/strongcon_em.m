clear all; 
AAPL1 = readmatrix('AAPL.csv');
AAPL = AAPL1(:,2); % read data
Xzero = AAPL(1);

randn('state', 100)
mu = 0.1244; sigma = 0.1038; % problem parameters
T = 1; N = 2^11; dt = T/N;
M = 10000; k = 5; % number of sample paths and test number of intervals

Xerr = zeros(M,k); Xplot = zeros(50,64); Xtrue1 = zeros(1,N);
for s = 1:M
    dW = sqrt(dt)* randn(1,N); % Brownian increments
    W = cumsum(dW); % discrete Brownian path
    Xtrue = Xzero*exp((mu-0.5*sigma^2)+sigma*W(end));
    for p = 1:k 
        R = 2^(p-1); Dt = R*dt; L = N/R; % L Euler steps of size Dt = R*dt 
        Xtemp = Xzero; 
        for j = 1:L 
            Winc = sum(dW(R*(j-1)+1:R*j)); 
            %Xtemp = Xtemp + Dt*mu*Xtemp + sigma*Xtemp*Winc; % EM method
            Xtemp = Xtemp + mu*Xtemp*Dt + sigma*Xtemp.*Winc+ 0.5*sigma^2*Xtemp.*(Winc.^2 - Dt); % Milstein method
            if (s<51) && (p==k)
                Xplot(s,j)= Xtemp;
            end
        end
        Xerr(s,p) = abs(Xtemp - Xtrue); % store the error at t = 1
    end
end

dt1 = 1/N; % interval of Xtrue
for i = 1:N
    Xtrue1(i) = Xzero*exp((mu-0.5*sigma^2)*dt1*i+sigma*W(i)); 
end



for i = 1:50
    plot([0:Dt:T], [Xzero, Xplot(i,:)],'r-'); hold on
end
h1 = plot([0:Dt: T], [Xzero, mean(Xplot)],'b-','LineWidth',3, 'DisplayName','Xmean'); hold on
h2 = plot([0:dt1:T], [Xzero , Xtrue1], 'g-' ,'LineWidth',2,'DisplayName','Xtrue'); hold off
legend([h1 h2]); 

mean(Xerr)

Dtvals = dt*(2.^([0:4]));
loglog(Dtvals,mean(Xerr),'b*-'), hold on
loglog(Dtvals,(Dtvals),'r--'), hold off % reference slope of 1/2
axis([1e-3 1e-1 1e-4 1])
xlabel('\Delta t'), ylabel('Sample average of | X(T) - X_L |')
% title('strongcon_em.m','FontSize',10)
title('milstein1','FontSize',10)
%%%% Least squares fit of error = C * Dt^q %%%%
A = [ones(5,1), log(Dtvals)']; rhs = log(mean(Xerr)');
sol = A\rhs; q = sol(2)
resid = norm(A*sol - rhs)
