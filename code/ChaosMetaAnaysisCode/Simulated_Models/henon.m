function [x,y,a]=henon(N,level,a,b)
% From Chaos Decision Tree Algorithm
% DOI: doi.org/10.6084/m9.figshare.7476362.v7

% Simulate the Henon map described in Henon (1976), "A two-dimensional
% mapping with a strange attractor"

%%% Inputs
% N - number of time-points to simulate
% level - level of white noise to add to the final signals relative to the standard deviation of those signals 
% each component of the Henon map 
% a, b - parameters of the map. 
% For periodic dynamics from the paper,
% set a=0.95, b=0.3

%
% Henon M (1976):A two-dimensional mapping with a strange attractor. 
% Communications in Mathematical Physics 50: 69-77

% Random initial conditions
x0=0.1*randn;
y0=0.1*randn;
x(1,1)=1-a*x0^2+b*y0;
y(1,1)=b*x0;

% Simulate
for i=2:N
    x(i,1)=1-a*x(i-1,1)^2+y(i-1,1);
    y(i,1)=b*x(i-1,1);
end

a=x+y;

% Add white noise
x=x+randn(N,1)*level*std(x);
y=y+randn(N,1)*level*std(y);
a=a+randn(N,1)*level*std(a);