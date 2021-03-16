function [x,y,z,R] =NPZ_model(N, level,regime,process_noise)

% Simulate an NPZ system described in Edwards et al. 1999,
% "The stability of an NPZ model subject to realistic levels of vertical
% mixing" 

%%% Inputs
% N - number of time-points to simulate
% level - level of lognormal noise to add to the final signal
% regime - a string to determine whether the script simulates the periodic
% ('periodic') or chaotic ('chaotic') dynamics
% process_noise - level of process noise to be added to a parameter value

%%% Outputs
% Time series

% Sample input: NPZ_model(100,0.1,'periodic',1)


if strmatch(regime, 'periodic')
    vm=2;
    ks=0.1;
    k=0.06;
    Rm=0.5;
    lam=0.2;
    gamma=0.3;
    m=0.1;
    G=0.2;
    sample_spacing=400;
    amp=0;
    depth=35;
elseif strmatch(regime, 'chaotic')
    vm=2;
    ks=0.1;
    k=0.06;
    Rm=4;
    lam=0.3;
    gamma=0.7;
    m=0.1;
    G=0.2;
    amp=1;
    sample_spacing=300;
    depth=0;
end


% random initial conditions
x(1)=3+0.2*randn;
y(1)=3;
z(1)=2;
noise_scale=0.05;

% time points, with an integration step of 0.01
h=.01;   %step size
t=0: h:  (sample_spacing*(N-1)*h+50000*h);
%t=0: h:  N;

% ordinary differential equations
f=@(t,x,y,z,R)  -vm*(1-amp*sin(2*pi/365*t))*x*y/(ks+x)*exp(-k*depth)+gamma*R*z*(1-exp(-lam*y))+m*y+G*z;
g=@(t,x,y,z,R) vm*(1-amp*sin(2*pi/365*t))*x*y/(ks+x)*exp(-k*depth)-R*z*(1-exp(-lam*y))-m*y;
p=@(t,x,y,z,R) (1-gamma)*R*z*(1-exp(-lam*y))-G*z;

% Simulate with RK4
for i=1:(length(t)-1)
    R(i)=Rm+process_noise*noise_scale*randn;
    
   k1=f(t(i),x(i),y(i),z(i),R(i));
    l1=g(t(i),x(i),y(i),z(i),R(i));
    m1=p(t(i),x(i),y(i),z(i),R(i));
    
      k2=f(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)),(R(i)+h/2));     
      l2=g(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)),(R(i)+h/2));
      m2=p(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)),(R(i)+h/2));
      
      k3=f(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)),(R(i)+h/2));
      l3=g(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)),(R(i)+h/2));
      m3=p(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)),(R(i)+h/2));
      
      k4=f(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h),(R(i)+h));
      l4=g(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h),(R(i)+h));
      m4=p(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h),(R(i)+h));
      
      x(i+1) = x(i) + h*(k1 +2*k2  +2*k3   +k4)/6; 
      y(i+1) = y(i) + h*(l1  +2*l2   +2*l3    +l4)/6;
      z(i+1) = z(i) + h*(m1+2*m2 +2*m3  +m4)/6;

end

% discard initial settling period
x=x(50001:sample_spacing:end);
y=y(50001:sample_spacing:end);
z=z(50001:sample_spacing:end);

% add  noise
l=length(x);
x=x'.*exp(normrnd(-(level*std(x))^2/2,level*std(x),l,1));
y=y'.*exp(normrnd(-(level*std(y))^2/2,level*std(y),l,1));
z=z'.*exp(normrnd(-(level*std(z))^2/2,level*std(z),l,1));