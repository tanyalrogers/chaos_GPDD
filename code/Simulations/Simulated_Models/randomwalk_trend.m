function y = randomwalk_trend(N,level,b)

% From Chaos Decision Tree Algorithm
% DOI: doi.org/10.6084/m9.figshare.7476362.v7

% Simulate a random walk 

%%% Inputs
% N - number of time-points to simulate
% level - level of observation noise 
% b - slope of the linear trend
% multiplied randn by 0.01 to generate data

%%% Outputs
% Time series

y=zeros(1,N);
y(1)=randn;
for i = 2:N
    y(i) = y(i-1)+b+randn;
end
y=y';
y=y+randn(N,1)*level*std(detrend(y));