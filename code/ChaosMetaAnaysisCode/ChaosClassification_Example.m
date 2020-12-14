%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script gives an example of generating data 
% and getting classification of 'chaotic' vs. 'not chaotic'
% for the methods in Supplementary Text B.3-B.6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; 

%% Declare variables
N = 250;  % length of time series
transient = 500; % number of time points to remove transients
level = 0.1;  % level of observation noise relative to standard deviation
E = 4; % embedding dimension (in paper, this was found using methods in Supplementary Text A.3)
tau = 1; % time lag (in paper, this was found using methods in Supplementary Text A.3)

%% Generate time series 
LongTimeSeries = logistic(N+transient,level,3.9); % call any function in 'Simulated_Models'
TimeSeries = LongTimeSeries(transient+1:end); % remove transients

%% Classify 
Classifications=ChaosClassification_MethodsB3toB6(TimeSeries,E,tau)