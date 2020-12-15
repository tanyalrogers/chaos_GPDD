function [x,y,a]=ikeda(n,level,mu)
% From Chaos Decision Tree Algorithm
% DOI: doi.org/10.6084/m9.figshare.7476362.v7

% Simulate the Ikeda map described in Hammel et al (1985), "Global
% dynamical behavior of the optical field in a ring cavity", based on
% earlier work in Ikeda (1979), "Multiple-valued stationary state and its
% instability of the transmitted light by a ring cavity system"

%%% Inputs
% N - number of time-points to simulate
% level - level of white noise to add to the final signals relative to the standard deviation of those signals 
% each component of the Ikeda map 
% mu - parameter of the map (set mu=0.9 for chaos as in paper)

%%% Outputs
% Time series

% Reference:
%
% Ikeda K (1979): Multiple-valued stationary state and its instability of the
% transmitted light by a ring cavity system. Optics Communications 30: 257

% random initial conditions
x0=0.1*randn;
y0=0.1*randn;
t=0.4-6/(1+x0^2+y0^2);
x(1,1)=1+mu*(x0*cos(t)-y0*sin(t));
y(1,1)=mu*(x0*sin(t)+y0*cos(t));

% Simulate
for i=2:n
    t=0.4-6/(1+x(i-1,1)^2+y(i-1,1)^2);
    x(i,1)=1+mu*(x(i-1,1)*cos(t)-y(i-1,1)*sin(t));
    y(i,1)=mu*(x(i-1,1)*sin(t)+y(i-1,1)*cos(t));
end

a=x+y;

% Add normal white noise
x=x+randn(n,1)*level*std(x);
y=y+randn(n,1)*level*std(y);
a=a+randn(n,1)*level*std(a);
