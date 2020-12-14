function [phi,beats]=poincare_oscillator(N,level,b,tau)
% From Chaos Decision Tree Algorithm
% DOI: doi.org/10.6084/m9.figshare.7476362.v7

% The following is slightly modified from the script originally
% written by Leon Glass for Nonlinear Dynamics in Physiology and Medicine
% (2003)

% This script simulates a periodically stimulated 1D Poincare oscillator

%%% Inputs
% N - number of time-points to simulate
% level - the amplitude of white noise to add to the final signals relative to the standard deviation of those signals 
% b, tau - parameters of the oscillator. b is the stimulation strength, and 
% tau is the period of the stimulation. 
% For chaotic dynamics, set b=1.13 and tau=0.65. 

%%% Outputs
% phi is a time-series of successive phases of the oscillator
% beats lists the number of beats between successive stimuli

%   Copyright Leon Glass 2003
%   Centre for Nonlinear Dynamics in Physiology and Medicine

phi=zeros(1,N);
phi(1)=rand;

for i=2:N+100
    angle= 2*pi*phi(i-1);
    rprime=sqrt(1+b^2+2*b*cos(angle));
    argument=(cos(angle)+b)/rprime;
    phi(i)=acos(argument)/(2*pi);
    if phi(i-1) > 0.5
        phi(i)=1-phi(i);
    end;
    phi(i)=phi(i)+tau;
    beats(i)=phi(i)-rem(phi(i),1);
    phi(i)=rem(phi(i),1);
end;

phi=phi(101:end)';
beats=beats(101:end)';

% Add white noise
phi=phi+randn(N,1)*level*std(phi);
beats=beats+randn(N,1)*level*std(beats);
