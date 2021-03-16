function [x,y] = SeasonalPredatorPrey(N, level,regime,process_noise)

% Simulate a seasonal predator prey system described in Turchin & Hanski (1997),
% "An empirically based model for latitudinal gradients in vole population
% dynamcis"

%%% Inputs
% N - number of time-points to simulate
% level - level of lognormal noise to add to the final signal
% regime - a string to determine whether the script simulates the periodic
% ('periodic') or chaotic ('chaotic') dynamics
% process_noise - level of process noise to be added to a parameter value

%%% Outputs
% Time series

% Sample input: SeasonalPredatorPrey(100,0.1,'periodic',1)

if strmatch(regime, 'periodic')
    r_m = 7;
    e=2.2;
    s = 1.5; 
    G=2;
    H=0.13;
    a=7.5;
    d=0.06;

elseif strmatch(regime, 'chaotic')
    r_m = 6;
    e=1;
    s = 1.25; 
    G=0;
    H=0.08;
    a=10;
    d=0.04;

end


% random initial conditions
x(1)=.10*rand;
y(1)=.10*rand;
noise_scale=0.5;

% time points, with an integration step of 0.01
h=.02;   %step size
t=0: h:  (10*(N-1)*h+  150000*h);

% ordinary differential equations
f=@(t,x,y,r) r*(1-e*sin(2*pi*t))*x-r*x^2-G*x^2/(x^2+H^2)-a*x*y/(x+d) ;
g=@(t,x,y) s*(1-e*sin(2*pi*t))*y-s*y^2/x;

r=r_m*ones(length(t)-1,1);

% Simulate with RK4
for i=1:(length(t)-1)
      r(i)=r_m+process_noise*noise_scale*randn;
    
      k1=f(t(i),x(i),y(i),r(i));
      l1=g(t(i),x(i),y(i));
    
      k2=f(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(r(i)+(0.5*h)));     
      l2=g(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)));
      
      k3=f(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(r(i)+(0.5*h)));
      l3=g(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)));
      
      k4=f(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(r(i)*h));
      l4=g(t(i)+h,(x(i)+k3*h),(y(i)+l3*h));

      x(i+1) = x(i) + h*(k1 +2*k2  +2*k3   +k4)/6; %final equations
      y(i+1) = y(i) + h*(l1  +2*l2   +2*l3    +l4)/6;

end

% discard initial settling period
x=x(150001:10:end);
y=y(150001:10:end);

% add noise
l=length(x);
x=x'.*exp(normrnd(-(level*std(x))^2/2,level*std(x),l,1));
y=y'.*exp(normrnd(-(level*std(y))^2/2,level*std(y),l,1));