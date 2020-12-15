function output = ChaosClassification_MethodsB3toB6(TimeSeries,E,tau)
% This script gets chaos classifications using the methods from 
% Supplementary Text B.3 - B.6
%
%%% INPUTS
% TimeSeries - time series of population dynamics
% E - the embedding dimension to be used in RQA
% tau - the time laag to be used in RQA
%
%%% OUTPUTS
% RecurrenceClass - Classification of 'chaotic' or 'not chaotic' from RQA
% PermutationClass - Classification of 'chaotic' or 'not chaotic' from permutation entropy
% VisibilityClass - Classification of 'chaotic' or 'not chaotic' from the horizontal visibility algorithm
% DecisionTreeClass - Classification of 'chaotic' or 'not chaotic' from the chaos decision tree
% Optional outputs: 
% DET - determinism from RQA
% L -  average length of diagonal line from RQA
% ENT - entropy from RQA
% PermutationEntropy - permutation entropy 
% MSE - deviation from expected degree distribution from the horizontal visibility algorithm
% VisibilityEntropy - entropy of the degree distribution from the horizontal visibility algorithm

n=length(TimeSeries);
 
%% RECURRENCE QUANTIFICATION ANALYSIS (Supplementary Text B.3)

% Source: 
% Cross Recurrence Plot (CRP) Toolbox (Ver 5.22 (R32.4)) (Marwan, 2020)
% THE CODE IN THIS SECTION NEEDS TO BE REQUESTED SEPARATELY
% GO TO: https://tocsy.pik-potsdam.de/CRPtoolbox/

% Call the RQA function
rqa_out=crqa(TimeSeries,E,tau,0.2,'euclidean','nonormalize','silent'); 

% Extract the most useful RQA measures 
DET=rqa_out(2);     % Determinism   
L =rqa_out(3);        % Mean length of diagonal lines
ENT =rqa_out(5);    %Shannon entropy of diagonal line lengths

%% PERMUTATION ENTROPY (Supplementary Text B.4)

% Source: 
% Toolboxes for Complex Systems (TOCSY) petropy function (Muller, 2019)
% The petropy function can be found at: https://tocsy.pik-potsdam.de/petropy.php

% Define the word length as 3 for everything because PE is badly biased for
% short time series 
wordLength=3; 

% Call the permutation entropy function 
PermutationEntropy=petropy(TimeSeries,wordLength,1);

%% HORIZONTAL VISIBILITY ALGORITHM (Supplementary Text B.5)

% Source: 
% Fast Horizontal Visibility Graph (HVG) for MATLAB file exchange (Iacobello, G., 2020)
% The fast_HVG function can be found at: 
% https://www.mathworks.com/matlabcentral/fileexchange/72889-fast-horizontal-visibility-graph-hvg-for-matlab

% Call fast horizontal visibility graph with default shift value
VisGraph=fast_HVG(TimeSeries,(1:n)',0);

% Get the degree of each node in visibility graph (i.e. point in time series)
Degrees=sum((VisGraph));

% Get the distribution of degrees 
maxDegree=20;
for k=2:maxDegree
    Pk(k-1)=sum(Degrees==k);
end

% Find indices where degree distribution is not zero
nonzeros=find(Pk~=0);

% Get the expected degree distribution for uncorrelated noise
expect=(1/3)*(2/3).^(0:maxDegree-2);

% Calculate the deviation from the expected degree distribution (scaled mean squared error)
% and entropy of degree distribution (scaled by log(4) to make it similar to LE) 
MSE=sqrt(mean((full(Pk(nonzeros))/(n-1)-expect(nonzeros)).^2))./std(expect(nonzeros));
VisibilityEntropy = full(- sum((Pk(nonzeros)/n).*log(Pk(nonzeros)/n)))-log(4);

%% DECISION TREE (Supplementary Text B.6)

% Source: 
% Chaos Decision Tree algorithm (Toker 2019)
% The chaos function can be found at: https://figshare.com/s/80891dfb34c6ee9c8b34

% Initialize all decision tree inputs
cutoff=[];
stationarity_test='bvr';
denoising_algorithm='schreiber';
gaussian_transform=0;
surrogate_algorithm='aaft_cpp';
downsampling_method='downsample';
sigma=0.0;

% Call chaos decision tree function with default values
DecisionTree=chaos(TimeSeries,cutoff,stationarity_test,denoising_algorithm,gaussian_transform,surrogate_algorithm,downsampling_method,sigma);


%% OUTPUT RESULTS

% Initialize all of the classification outputs to 'not chaotic'
RecurrenceClassification='not chaotic';
PermutationClassification='not chaotic';
VisibilityClassification='not chaotic';
DecisionTreeClassification='not chaotic';

% If the appropriate criteria are met, classify as chaotic
if (DET>0.45 && DET < .99) &&(L>1.9 && L < 5.3) && (ENT>0.39 && ENT < 2.3)
    RecurrenceClassification='chaotic';
end 

if (PermutationEntropy>1.06 && PermutationEntropy < 1.23) 
    PermutationClassification='chaotic';    
end 

if (VisibilityEntropy>0.32 && VisibilityEntropy < .46) && (MSE>0.21 && MSE < .48)
    VisibilityClassification='chaotic';
end 

if strcmp(DecisionTree.result,'chaotic')
    DecisionTreeClassification='chaotic';
end 

% Output all of the classifications and important measures (optional)
output.RecurrenceClass=RecurrenceClassification; 
output.PermutationClass=PermutationClassification;
output.VisibilityClass=VisibilityClassification; 
output.DecisionTreeClass=DecisionTreeClassification; 

% % OPTION: output the important measures (uncomment to include these in output) 
% output.DET=DET;
% output.L=L;
% output.ENT=ENT;
% output.PermutationEntropy=PermutationEntropy;
% output.MSE=MSE; 
% output.VisibilityEntropy=VisibilityEntropy; 
end 