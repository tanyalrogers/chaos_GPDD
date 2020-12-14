function y = cyclostationary_AR(N,tau,T)
% From Chaos Decision Tree Algorithm
% DOI: doi.org/10.6084/m9.figshare.7476362.v7

% Simulate the cyclostationary autoregressive process described in Timmer,
% 1998, "Power of surrogate data testing with respect to nonstationarity."

%%% Inputs
% N - number of time-points to simulate
% tau - relaxation time
% T - oscillation period

%%% Parameters
%tau=50;
%T=10;

%%% Outputs
% Time series

a1=2*cos(2*pi/T).*exp(-1/tau);
a2=-exp(-2/tau);
y=zeros(1,N);
y(1:2)=randn;
for i = 3:N
    y(i) = a1*y(i-1)+a2*y(i-2)+randn;
end