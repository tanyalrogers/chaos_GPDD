function [x,y]=predatorprey(n, level, regime)

% Simulate predator prey system from Zhang et al. (2018) " Complex Dynamics 
% on the Routes to Chaos in a Discrete Predator-Prey System with Crowley-Martin 
% Type Functional Response"

%%% Inputs
% n - number of time-points to simulate
% level - level of lognormal noise to add to the final signal relative to the standard deviation
% regime - a string to determine whether the script simulates the periodic
% ('periodic') or chaotic ('chaotic') dynamics

%%% Outputs
% Time series

% Sample input: predatorprey(100,0.1,'periodic')

% set up parameters for the specified regime
if strmatch(regime, 'periodic')
    a = 2;
    b=2;
    c=2; 
    d=1.85;
    alpha = 0.1; 
    beta = 0.1; 
    tau=1.1;
elseif strmatch(regime, 'chaotic')
    a = 2;
    b=2;
    c=2; 
    d=1.85;
    alpha = 0.1; 
    beta = 0.1; 
    tau=1.27;
end

% Random initial conditions within reasonable range
x(1)=1+rand;
y(1)=0.2+0.4*rand;

% Simulate
for i=2:n
    x(i)=x(i-1)+tau*x(i-1)*(a-x(i-1)-(b*y(i-1))/((1+alpha*x(i-1))*(1+beta*y(i-1)))); 
    y(i)=y(i-1)+tau*y(i-1)*(-c+(d*x(i-1))/((1+alpha*x(i-1))*(1+beta*y(i-1)))); 
end

std_x=std(x);

% Add log-normal noise
scale=0.15+0.65*exp(-std_x);
x=x'.*lognrnd(0,scale*level*std_x,n,1);
y=y';



