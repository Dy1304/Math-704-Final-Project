clear all; clc;

AAPL1 = readmatrix('AAPL.csv'); % read data from AAPL.csv
AAPL = AAPL1(:,2);
Xzero = AAPL(1);
n = length(AAPL); Dt = 1/n; T=1;
plot([Dt:Dt:T], AAPL,'r-') % plot the opening price of AAPL in the stock market
xlabel('t')
ylabel('AAPL')
xlim([0 1])
ylim([110 190])

sum0 = 0;
for i = 1:n-1
    sum0 = sum0 + (AAPL(i+1)/AAPL(i)) - 1;
end
mu = sum0/T; % obtain the ME of mu

sum1 = 0;
for i = 1:n-1
    sum1 = sum1 + ((AAPL(i+1)/AAPL(i)) - 1 - mu*(1/n))^2;
end
sigma = sum1/T; % obtain the ME of sigma