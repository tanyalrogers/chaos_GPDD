function [x,y]=tinkerbellmap(n, level, regime)

% Simulate Tinkerbell map from Shaoliang et al. (2011) "Bifurcation and chaos in the Tinkerbell map"

%%% Inputs
% n - number of time-points to simulate
% level - level of lognormal noise to add to the final signal relative to the standard deviation
% regime - a string to determine whether the script simulates the periodic
% ('periodic') or chaotic ('chaotic') dynamics

%%% Outputs
% Time series

% Sample input: tinkerbellmap(100,0.1,'periodic')

% set up parameters for the specified regime
if strmatch(regime, 'periodic')
    a =0.9;
    b=-0.5;
    c=1.8; 
    d=0.5;
elseif strmatch(regime, 'chaotic')
    a =0.9;
    b=-0.5;
    c=2.15; 
    d=0.5;
end

% Random initial conditions within reasonable range
x(1)=-0.05*rand;
y(1)=-0.05*rand;

% Simulate
for i=2:n
    x(i)=x(i-1)^2-y(i-1)^2+a*x(i-1)+b*y(i-1); 
    y(i)=2*x(i-1)*y(i-1)+c*x(i-1)+d*y(i-1); 
end


% Add lognormal noise
x=(x'+1).*lognrnd(0,level*std(x),n,1);
y=y';