function [x,theta,y]=cubicmap(n, level, regime)
% Adapted from Chaos Decision Tree Algorithm
% DOI: doi.org/10.6084/m9.figshare.7476362.v7

% Simulate the periodically forced cubic map described in Venkatesan and 
% Lakshmanan (2001), "Interruption of torus doubling bifurcation and 
% genesis of strange nonchaotic attractors in a quasiperiodically forced 
% map: Mechanisms and their characterizations"

%%% Inputs
% N - number of time-points to simulate
% level - level of white noise to add to the final signals relative to the standard deviation of those signals 
% regime - a string to determine whether the script simulates the periodic
% ('periodic') or chaotic ('chaotic') dynamics

%%% Outputs
% Time series

% Sample input: cubicmap(100,0.1,'periodic')

w = (sqrt(5)-1)/2; % golden ratio

% Random initial conditions
x(1,1)=abs(randn);
theta(1,1)=abs(randn);

% set up parameters for the specified regime
if strmatch(regime, 'periodic')
    f = 0;
    Q = 0;
    A = 1;
elseif strmatch(regime, 'chaotic')
    f=-0.8;
    Q = 0;
    A = 1.5;
elseif strmatch(regime,'HH')
    % SNA through the Heagy-Hammel route
    A = 1.88697;
    Q = 0;
    f = 0.7;
elseif strmatch(regime, 'S3')
    % SNA through type-3 intermittency
    A = 2.14;
    Q = 0;
    f = 0.35;
elseif strmatch(regime, '2T')
    A=1.1;
    Q=0;
    f=-0.18;
end


% Simulate
for i=2:n+100
    
    x(i,1)=Q + f*cos(2*pi*theta(i-1,1)) - A*x(i-1,1) + x(i-1,1)^3;

    theta(i,1)=mod((theta(i-1,1)+w),1);
end
x=x(101:end);
theta=theta(101:end);
% Take linear combination of x and theta
y = x./(6)+theta./10; 

% Add white noise
x=x+randn(n,1)*level*std(x);
theta=theta+randn(n,1)*level*std(theta);
y=y+randn(n,1)*level*std(y);
