%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script applies RQA, RE, HVG, and CDT methods to all simulated data 
% in the test dataset and the validation dataset
%
% NOTE: This script takes >20 hours to run in full
%
% Results from running the TEST DATA section are in 
% "sims_test_results_othermethods.csv" in the GitHub data folder
%
% Results from running the VALIDATION DATA section are in 
% "sims_validation_results_othermethods.csv" in the GitHub data folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; 

%% TEST DATA 

% Import simulated data with corresponding E and tau values
test_data=readtable('simulation_dataset_test.csv');
Es=readtable('sims_test_E.csv');
taus=readtable('sims_test_tau.csv');

% Preallocate output for speed
ID=cell(length(unique(test_data.ID))*100,1);
SimNumber=cell(length(unique(test_data.ID))*100,1);
RQAclass=cell(length(unique(test_data.ID))*100,1);
PEclass=cell(length(unique(test_data.ID))*100,1);
HVAclass=cell(length(unique(test_data.ID))*100,1);
DTclass=cell(length(unique(test_data.ID))*100,1);

% Apply RQA, PE, HVG, and CDT methods to simulated data. 
%NOTE: This takes a long time to run!
count=1; % initialize loop count

for j = 1:100
    for i = 1:length(unique(test_data.ID)) 
        
        % Get time series, E, and tau
        TimeSeries= table2array(test_data(test_data.ID==i,join(['Sim', sprintf("%d",j)],"_"))); 
        E = table2array(Es(Es.ID==i,join(['Sim', sprintf("%d",j)],"_")));
        tau = table2array(taus(taus.ID==i,join(['Sim', sprintf("%d",j)],"_")));
        
        % Apply the classification function
        out=ChaosClassification_MethodsB3toB6(TimeSeries,E,tau);
        
        % Store results for each method
        ID{count} = i;
        SimNumber{count}=join(['Sim', sprintf("%d",j)],".");
        RQAclass{count}=out.RecurrenceClass;
        PEclass{count}=out.PermutationClass;
        HVAclass{count}=out.VisibilityClass;
        DTclass{count}=out.DecisionTreeClass;
        
        % Iterate count
        count=count+1;
    end 
end 

% Merge results into a table and write to csv
output = [cell2table(ID) cell2table(SimNumber) cell2table(RQAclass) cell2table(PEclass) cell2table(HVAclass) cell2table(DTclass)];
%writetable(output, 'sims_test_results_othermethods.csv')  % NOTE: "sims_test_results_othermethods.csv" is already in GitHub data folder


%% VALIDATION DATA 
clear all; close all; 

% Import simulated data with corresponding E and tau values
validation_data=readtable('simulation_dataset_validation.csv');
Es=readtable('sims_validation_E.csv');
taus=readtable('sims_validation_tau.csv');

% Preallocate output for speed
ID=cell(length(unique(validation_data.ID))*100,1);
SimNumber=cell(length(unique(validation_data.ID))*100,1);
RQAclass=cell(length(unique(validation_data.ID))*100,1);
PEclass=cell(length(unique(validation_data.ID))*100,1);
HVAclass=cell(length(unique(validation_data.ID))*100,1);
DTclass=cell(length(unique(validation_data.ID))*100,1);

% Apply RQA, PE, HVG, and CDT methods to simulated data. NOTE: This takes a long time to run!
count=1; % initialize loop count
 
for j = 1:100
    for i = 1:length(unique(validation_data.ID))   
        % Get time series, E, and tau
        TimeSeries= table2array(validation_data(validation_data.ID==i,join(['Sim', sprintf("%d",j)],"_"))); 
        E = table2array(Es(Es.ID==i,join(['Sim', sprintf("%d",j)],"_")));
        tau = table2array(taus(taus.ID==i,join(['Sim', sprintf("%d",j)],"_")));
        
        % Apply the classification function
        out=ChaosClassification_MethodsB3toB6(TimeSeries,E,tau);
        
        % Store results for each method
        ID{count} = i;
        SimNumber{count}=join(['Sim', sprintf("%d",j)],".");
        RQAclass{count}=out.RecurrenceClass;
        PEclass{count}=out.PermutationClass;
        HVAclass{count}=out.VisibilityClass;
        DTclass{count}=out.DecisionTreeClass;
        
        % Iterate count
        count=count+1;
    end 
end 

% Merge results into table and write to a csv file
output = [cell2table(ID) cell2table(SimNumber) cell2table(RQAclass) cell2table(PEclass) cell2table(HVAclass) cell2table(DTclass)];
%writetable(output, 'sims_validation_results_othermethods.csv')  % NOTE: "sims_validation_results_othermethods.csv" is already in GitHub data folder



%% NOISY DATA 
clear all; close all; 

% Import simulated data with corresponding E and tau values
noise_data=readtable('simulation_dataset_noise_test.csv');
Es=readtable('sims_noise_E.csv');
taus=readtable('sims_noise_tau.csv');

% Preallocate output for speed
ID=cell(length(unique(noise_data.ID))*100,1);
SimNumber=cell(length(unique(noise_data.ID))*100,1);
RQAclass=cell(length(unique(noise_data.ID))*100,1);
PEclass=cell(length(unique(noise_data.ID))*100,1);
HVAclass=cell(length(unique(noise_data.ID))*100,1);
DTclass=cell(length(unique(noise_data.ID))*100,1);

% Apply RQA, PE, HVG, and CDT methods to simulated data. 
count=1; % initialize loop count
 
for j = 1:100
    for i = 1:length(unique(noise_data.ID))   
        % Get time series, E, and tau
        TimeSeries= table2array(noise_data(noise_data.ID==i,join(['Sim', sprintf("%d",j)],"_"))); 
        TimeSeries=(TimeSeries-mean(TimeSeries))/std(TimeSeries);
        E = table2array(Es(Es.ID==i,join(['Sim', sprintf("%d",j)],"_")));
        tau = table2array(taus(taus.ID==i,join(['Sim', sprintf("%d",j)],"_")));
        
        % Apply the classification function
        out=ChaosClassification_MethodsB3toB6(TimeSeries,E,tau);
        
        % Store results for each method
        ID{count} = i;
        SimNumber{count}=join(['Sim', sprintf("%d",j)],".");
        RQAclass{count}=out.RecurrenceClass;
        PEclass{count}=out.PermutationClass;
        HVAclass{count}=out.VisibilityClass;
        DTclass{count}=out.DecisionTreeClass;
        
        % Iterate count
        count=count+1;
    end 
end 

% Merge results into table and write to a csv file
output = [cell2table(ID) cell2table(SimNumber) cell2table(RQAclass) cell2table(PEclass) cell2table(HVAclass) cell2table(DTclass)];
%writetable(output, 'sims_noise_results_othermethods.csv')  % NOTE: "sims_noise_results_othermethods.csv" is already in GitHub data folder