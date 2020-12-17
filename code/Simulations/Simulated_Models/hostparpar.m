function [x,y]=hostparpar(n, level, regime)

% Simulate host-parasitoid-parasitoid system from Yu et al. (2009) " Dynamic complexities in
% a parasitoid-host-parasitoid ecological model"

%%% Inputs
% n - number of time-points to simulate
% level - level of lognormal noise to add to the final signal relative to the standard deviation
% regime - a string to determine whether the script simulates the periodic
% ('periodic') or chaotic ('chaotic') dynamics

%%% Outputs
% Time series

% Sample input: hostparpar(100,0.1,'periodic')

% set up parameters for the specified regime
if strmatch(regime, 'periodic')
    beta=1.4; 
    K=20;
    a=0.9;
    b=0.8; 
    mm=0.7;  
    nn=0.4;
    r=2.5; 
elseif strmatch(regime, 'chaotic')
    beta=1.4; 
    K=20;
    a=0.9;
    b=0.8; 
    mm=0.7;  
    nn=0.4;
    r=3.4; 
end

% Random initial conditions
x(1)=rand;
y(1)=rand;
z(1)=rand;

% Simulate
for t=1:n-1
     x(t+1)= (x(t).*exp(r.*(1-x(t)/K)-a*y(t).^(-mm+1)-b*beta.*z(t).^(-nn+1)));
     y(t+1)= (x(t).*(1-exp(-a*y(t).^(-mm+1)-b*beta.*z(t).^(-nn+1))).*(a*y(t).^(-mm+1)./(a*y(t).^(-mm+1)+b*beta.*z(t).^(-nn+1))));
     z(t+1)= (x(t).*(1-exp(-a*y(t).^(-mm+1)-b*beta.*z(t).^(-nn+1))).*(b*beta*z(t).^(-nn+1)./(a*y(t).^(-mm+1)+b*beta.*z(t).^(-nn+1))));
end

std_x=std(x);

% Add lognormal noise
scale=0.15+0.65*exp(-std_x); % lognormal noise can get really big, so we scale
x=x'.*lognrnd(0,scale*level*std_x,n,1);
y=y';
z=z';
