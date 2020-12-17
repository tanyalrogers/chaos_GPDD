function [x,y]=threespeciescompetition(n, level, regime)

% Simulate competition dynamics from Ali et al. (2019) "Bifurcation analysis and chaos control
% in discrete-timesystem of three competing species"

%%% Inputs
% n - number of time-points to simulate
% level - level of lognormal noise to add to the final signal relative to the standard deviation
% regime - a string to determine whether the script simulates the periodic
% ('periodic') or chaotic ('chaotic') dynamics

%%% Outputs
% Time series

% Sample input: threespeciescompetition(100,0.1,'periodic')

% set up parameters for the specified regime
if strmatch(regime, 'periodic')
    a=0.65;
    b=0.6; 
    r=2.6; 
elseif strmatch(regime, 'chaotic')
    a=0.65;
    b=0.6; 
    r=3; 
end

% Random initial conditions
x(1)=rand;
y(1)=rand;
z(1)=rand;

% Simulate
for t=1:n-1
     x(t+1)= (x(t).*exp(r.*(1-x(t)-a*y(t)-b*z(t))));
     y(t+1)= (y(t).*exp(r.*(1-y(t)-a*z(t)-b*x(t))));
     z(t+1)= (z(t).*exp(r.*(1-z(t)-a*x(t)-b*y(t))));
end

std_x=std(x);

% Add log-normal noise
scale = 0.15+0.65*exp(-std(x));
x=x'.*lognrnd(0,scale*level*std_x,n,1);
y=y';
z=z';