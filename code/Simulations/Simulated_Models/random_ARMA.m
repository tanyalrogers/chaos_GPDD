function y = random_ARMA(N, level, p,theta)

% random ARMA(1) model with Gaussian-distributed error

y(1)=randn;
err=randn(1,N);
for i = 2:N
    y(i) = p*y(i-1)+err(i)+(1+theta)*err(i-1);
end
y=y';
y=y+randn(N,1)*level*std(y);