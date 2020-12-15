function y = AR1(N,level)
% Simulate an autoregressive model

%%% Inputs
% N - number of time-points to simulate
% level - noise level

%%% Outputs
% time series

y=zeros(1,N);
y(1)=randn;
for i = 2:N
    y(i) = 8-0.8*y(i-1)+randn;
end
y=y';

% Add white observation noise
y=y+randn(N,1)*level*std(y);