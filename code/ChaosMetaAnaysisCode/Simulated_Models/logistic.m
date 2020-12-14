function y=logistic(N,level,r,y0)
% Adapted from Chaos Decision Tree Algorithm
% DOI: doi.org/10.6084/m9.figshare.7476362.v7

% Simulate the logistic map described in May (1976), "Simple mathematical
% models with very complicated dynamics"

%%% Inputs
% N - number of time-points to simulate
% level - level of white noise to add to the final signal relative to the standard deviation of the signals 
% r - growth rate parameter
% y0 - initial value (optional)

%%% Outputs
% Time series

%%% Other information
% r = 3.55 for 8-cycle
% r = 3.828427 for 3-cycle
% r=3.9 for chaos
 
% random initial condition
if nargin<4
y0=rand;
end
y(1,1)=r*y0*(1-y0);

% Simulate
for i=2:N
    y(i,1)=r*y(i-1,1)*(1-y(i-1,1));
end

y=y+randn(N,1)*level*std(y);

