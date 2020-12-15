function y = AR2(N,level)
% Simulate an autoregressive model

%%% Inputs
% N - number of time-points to simulate
% level - noise level

%%% Outputs
% time series

y=zeros(1,N);
y(1)=rand;
y(2)=rand;
for i = 3:N
    y(i) = 0.9*y(i-1)-0.1*y(i-2)+randn;
end
y=y';

% Add white observation noise
y=y+randn(N,1)*level*std(y);