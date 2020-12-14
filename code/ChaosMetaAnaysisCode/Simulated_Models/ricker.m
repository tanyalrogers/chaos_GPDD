function y=ricker(N,level,r,y0)

% Simulate Ricker model

%%% Inputs
% N - number of time-points to simulate
% level - level of white noise to add to the final signal relative to the standard deviation of the signals 
% r - growth rate parameter
% y0 - initial value (optional)

%%% Outputs
% Time series

%%% Other information
% r = 2.2 for 2-cycle
% r=3.4 for chaos

% random initial condition
if nargin<4
y0=rand;
end
y(1,1)=y0*exp(r*(1-y0));

% Simulate
for i=2:N
    y(i,1)=y(i-1,1)*exp(r*(1-y(i-1,1)));
end

y=y+randn(N,1)*level*std(y);