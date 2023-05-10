T = 1; N = 200; dt=T/N; 

dW = sqrt(dt)*randn(1,N); % set each increment
W = cumsum(dW); % cumulative sum

plot([0:dt:T], [0,W], 'r--')
xlabel('t')
ylabel('W(t)')
