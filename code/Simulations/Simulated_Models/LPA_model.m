function [x,y,z]=LPA_model(n, level, regime,process_noise)

% Simulate LPA model from Constantino et al. (1997) 
% "Chaotic Dynamics in an Insect Population"

%%% Inputs
% n - number of time-points to simulate
% level - level of lognormal noise to add to the final signal
% regime - a string to determine whether the script simulates the periodic
% ('periodic') or chaotic ('chaotic') dynamics
% process_noise - level of process noise to be added to the dynamics

%%% Outputs
% Time series

% Sample input: LPA_model(100,0.1,'periodic',1)

% set up parameters for the specified regime
if strmatch(regime, 'periodic')
    mu_l= 0.2055;
    mu_a=0.96;
    c_el=0.01209;
    c_ea=0.01155;
    c_pa=1;
    b=6.598;
elseif strmatch(regime, 'chaotic')
    mu_l= 0.2055;
    mu_a=0.96;
    c_el=0.01209;
    c_ea=0.01155;
    c_pa=0.35;
    b=6.598;
end

% Random initial conditions within reasonable range
x(1)=rand;
y(1)=rand;
z(1)=rand;
noise_scale=0.035;

% Simulate
for i=2:n+500
    x(i)=b*z(i-1)*exp(-c_el*x(i-1)-c_ea*z(i-1))*exp(normrnd(-(process_noise*noise_scale)^2/2,(process_noise*noise_scale)));
    y(i)=(1-mu_l)*x(i-1)*exp(normrnd(-(process_noise*noise_scale)^2/2,(process_noise*noise_scale)));
    z(i)=(y(i-1)*exp(-c_pa*z(i-1))+(1-mu_a)*z(i-1))*exp(normrnd(-(process_noise*noise_scale)^2/2,(process_noise*noise_scale))); 
end

x=x(501:end);
y=y(501:end);
z=z(501:end);

% Add log-normal noise
x=x'.*exp(normrnd(-(level)^2/2,level,n,1));
y=y'.*exp(normrnd(-(level)^2/2,level,n,1));
z=z'.*exp(normrnd(-(level)^2/2,level,n,1));
