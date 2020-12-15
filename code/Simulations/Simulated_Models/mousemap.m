function [x,y]=mousemap(n, level, regime)

% Simulate mouse (Gauss) map from Hilborn (2004) "Chaos and nonlinear dynamics:
% an introduction for scientists and engineers"

%%% Inputs
% n - number of time-points to simulate
% level - level of lognormal noise to add to the final signal relative to the standard deviation
% regime - a string to determine whether the script simulates the periodic
% ('periodic') or chaotic ('chaotic') dynamics

%%% Outputs
% Time series

% Sample input: mousemap(100,0.1,'periodic')

% set up parameters for the specified regime
if strmatch(regime, 'periodic')
    alpha = 6.2; 
    beta = 0.0; 
elseif strmatch(regime, 'chaotic')
    alpha = 6.2; 
    beta = -0.5; 
end

% Random initial conditions
x(1)=rand;

% Simulate
for i=2:n
    x(i)=exp(-alpha*x(i-1).^2)+beta; 
end

% Add lognormal noise
x=x'.*lognrnd(0,level*std(x),n,1);