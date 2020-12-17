%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script applies RQA, RE, HVG, and CDT methods to empirical data 
% from the GPDD
%
% Results from running this script are in 
% "gpdd_results_othermethods.csv" in the GitHub data folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; 

% Import simulated data with corresponding E and tau values
gpdd_data=readtable('gpdd_timeseries.csv');
meta_data=readtable('gpdd_ts_metadata.csv');
EandTau=readtable('gpdd_Etau_smap.csv');

% Get MainIDs
MainIDs=meta_data.MainID;

% Preallocate output for speed
MainID=cell(length(MainIDs),1);
RQA=cell(length(MainIDs),1);
PE=cell(length(MainIDs),1);

for i = 1:length(MainIDs) 
    % MainID for this iteration
    IDnum=MainIDs(i);
    
    % Get the time series, E, and tau
     indices= find(gpdd_data.MainID==IDnum);
     TimeSeries =str2double(gpdd_data.PopRescale(indices));
     E=EandTau.E(EandTau.MainID==IDnum);
     tau=EandTau.tau(EandTau.MainID==IDnum);
        
     % Apply the classification function
     out=ChaosClassification_MethodsB3toB6(TimeSeries,E,tau);
        
     % Store results for each method
     MainID{i} = MainIDs(i);
     RQA{i}=out.RecurrenceClass;
     PE{i}=out.PermutationClass;
end 

% Merge results into a table and write to csv
output = [cell2table(MainID) cell2table(RQA) cell2table(PE)];
%writetable(output, 'gpdd_results_othermethods.csv')  % NOTE: "gpdd_results_othermethods.csv" is already in GitHub data folder
