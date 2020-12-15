function y=sineWave(N,level,a,period)

% Simulate sine wave

%%% Inputs
% N - number of time-points to simulate
% level - level of white noise to add to the final signal relative to the standard deviation 
% a - amplitude of sine wave
% period - period of sine wave

%%% Outputs
% Time series

%%% Other information
% a = 1
% b = 12

% random initial condition
y(1,1)=a*sin(2*pi/period*1);

% Simulate
for i=2:N
    y(i,1)=a*sin(2*pi/period*i)+a;
end

% Add white observation noise
y=y+randn(N,1)*level*std(y);