function output = chaos(y,cutoff,stationarity_test,denoising_algorithm,...
    gaussian_transform,surrogate_algorithm,downsampling_method,sigma)

% INPUTS

% y % 
    % a time-series vector

% cutoff % 
    % The cutoff for the K-statistic, outputed by the 0-1 test for chaos, 
    % for classifying data as chaotic or periodic. If no cutoff is
    % provided, then a cutoff will be automatically chosen based on the
    % length of the time-series (Supplementary Figure 4)

% stationarity_test %
    % If you want to first test your data for stationarity before proceeding
    % to test for stochasticity or chaos, set stationarity_test to a string
    % indicating which stationarity test you wish to use:
    %
    % 'adf' for the Augmented Dickey-Fuller test
    % 'kpss' for the  Kwiatkowski–Phillips–Schmidt–Shin test 
    % 'lmc' for the Leybourne-McCabe test
    % 'vrat' for Lo and Mckinlay's standard variance ratio test
    % 'bvr' for Breitung's variance ratio test 
    %
    % Among these, we found that Breitung's variance ratio test was most 
    % robust to edge cases (Supplementary Table 1). If this variable is 
    % empty, then no stationarity test will be used (default)

% denoising_algorithm %
    % A string indicating which among three denoising algorithms to use: 
    %
    % 'schreiber' for Schreiber's de-noising algorithm (default)
    % 'ma' for a moving average filter
    % 'wavelet' for wavelet denoising using an empirical Bayesian method 
    %  with a Cauchy prior

% gaussian_transform %
    % A flag indicating whether or not to normality-transform the data,
    % using the Box-Cox method (1), the rank-based inverse normal
    % transformation method (2), or no transformation (0) before generating 
    % surrogate data:
    %
    % 0 - no normality transformation (default)
    % 1 - normality-transform using Box-Cox method
    % 2 - normality-transform using rank-based inverse normal
    % transformation
    
    
% surrogate_algorithm %
    % A string indicating which surrogate algorithm to use for the
    % surrogate-based stochasticity test. Options are:
    %
    % 'IAAFT2' - iterated amplitude adjusted Fourier transform algorithm
    % 'AAFT' - amplitude adjusted Fourier transform algorithm
    % 'FT' - Fourier transform algorithm
    % 'CPP' - cyclic phase permutation algorithm
    % 'CSS' - cycle shuffle surrogate algorithm
    % 'aaft_cpp' - a combination of both AAFT and CPP surrogates (default)

% downsampling_method %
    % A string indicating which downsampling method to use:
    %
    % 'downsample' - standard downsampling (default)
    % 'fouda' - take the local minima and maxima of the signal, as
    % recommended by Fouda and colleagues, 2013, "A Modified 0-1 Test for 
    % Chaos Detection in Oversampled Time Series Observations"
    
% sigima %
    % A parameter controlling the level of noise used to suppress
    % correlations in the modified 0-1 test. Default = 0.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default tests and parameters if none are specified

h=0;  

if nargin <4 || isempty(denoising_algorithm)
    denoising_algorithm='schreiber';
end

if nargin <4 || isempty(gaussian_transform)
    gaussian_transform=0;
end

if nargin <5 || isempty(surrogate_algorithm)
    surrogate_algorithm = 'aaft_cpp';
end

if nargin <7 || isempty(downsampling_method)
    downsampling_method='downsample';
end

if nargin <8 || isempty(sigma)
    sigma=0.5;
end

num_surr = 1000; % number of surrogates time-series to generate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Optional test for non-stationarity 
if nargin>1 && ~isempty(stationarity_test)
    if strcmp(stationarity_test, 'adf')
        h=~adftest(y);
    elseif strcmp(stationarity_test, 'lmc')
        h=lmctest(y);
    elseif strcmp(stationarity_test, 'kpss')
        h=kpsstest(y);
    elseif strcmp(stationarity_test, 'vrat')
        h=~vratiotest(y);
    elseif strcmp(stationarity_test, 'bvr')
        h=~bvr(y,0);
    end
end

if h==1
    output.result='nonstationary';
else    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 2: Test if data are deterministic
    
    if gaussian_transform == 1
        surr_y = boxcox(y+abs(min(y))+.0001); % add offset so all values are positive
    elseif gaussian_transform == 2
        surr_y = rank_inverse_gaussian_transform(y);
    else
        surr_y=y;
    end
    
    if size(surr_y,1)>size(surr_y,2)
        surr_y=surr_y';
    end
    
    % Test if data are stochastic and linear using permutation entropy
    % and AAFT surrogates
    %%
    if strcmp(surrogate_algorithm, 'aaft_cpp')
        try
            [surr, params] = surrogate(surr_y, num_surr, 'AAFT', 1, 1);
        catch
            [surr, params] = surrogate(zscore(surr_y), num_surr, 'AAFT', 1, 1);
        end
        sig=params.cutsig;
        for i = 1:num_surr
            surr_h1(i) = petropy(surr(i,:),8,1);
        end
        perm_h1 = petropy(sig,8,1);
        
        stochastic1 = perm_h1>=min(surr_h1)&&perm_h1<=max(surr_h1);
        
        % Stochastic *nonlinear* data might pass the prior test and be classified
        % as deterministic. To try to rule this out, do a second test of determinism using
        % cyclic phase permutation surrogates, which maintains individual cycles in
        % the data but shuffles them, thus breaking determinism
        try
            [surr, params] = surrogate(surr_y, num_surr, 'CPP', 1, 1);
        catch
            [surr, params] = surrogate(zscore(surr_y), num_surr, 'CPP', 1, 1);
        end
        sig=params.cutsig;
        for i = 1:num_surr
            x1 = surr(i,:);
            surr_h2(i) = petropy(x1,8,1);
        end
        perm_h2 = petropy(sig,8,1);
        
        stochastic2 = perm_h2>=min(surr_h2)&&perm_h2<=max(surr_h2);
        
        stochastic = stochastic1 || stochastic2;
    else
        try
            [surr, params] = surrogate(surr_y, 1000, surrogate_algorithm, 1, 1);
        catch
            [surr, params] = surrogate(zscore(surr_y), 1000, surrogate_algorithm, 1, 1);
        end
        sig=params.cutsig;
        for i = 1:1000
            surr_h(i) = petropy(surr(i,:),8,1);
        end
        perm_h = petropy(sig,8,1);
        
        stochastic = perm_h>=min(surr_h)&&perm_h<=max(surr_h);
    end
    
    if stochastic
        output.result = 'stochastic';
    else
     %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 3: De-noise the data
        if strcmp(denoising_algorithm, 'schreiber')
            y=noiserSchreiber(y);
            y=y(10:end-10); % cut off ends
        elseif strcmp(denoising_algorithm, 'ma') % moving average
            y=smooth(y);
        elseif strcmp(denoising_algorithm, 'wavelet')
            y=wdennoise(y); % requires Matlab 2018a or more recent
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4: Check for over-sampling. If data are oversampled,
        % downsample until no longer over-sampled
        
        % crude test of over-sampling:
        if (max(y)-min(y) )/mean(abs(diff(y))) >10
            if strcmp(downsampling_method,'fouda')
                y=minmaxsig(y); % requires Matlab 2018a or more recent
            elseif strcmp(downsampling_method, 'downsample')
                oversample_flag=1;
                while oversample_flag==1
                    y=downsample(y,2);
                    if (max(y)-min(y) )/mean(abs(diff(y))) < 10 ...
                            || length(y)<100
                        oversample_flag=0;
                    end
                end
            end
        end
        
             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 6: Test for presence of chaos using modified 0-1 test,
        % and degree of chaos using permutation entropy
        
        % set K-statistic cutoff based on time-series length, if no cutoff
        % is provided
        if nargin <2 || isempty(cutoff)
            load('cutoff_fit.mat')
            cutoff=f(length(y));
            if cutoff>.99
                cutoff=.99;
            end
        end
        
        % permutation entropy for degree of chaos
        output.permutation_entropy = petropy(y,5,1);
        
        % Normalize standard deviation of signal before applying 0-1 test
        norm_fac=.5/std(y);
        y=y.*norm_fac;
        
        % modified 0-1 test
        K=z1test(y,sigma);
        if K>cutoff
            output.result = 'chaotic';
        else
            output.result = 'periodic';
        end
        output.K = K;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code for calculating permutation entropy is written by
% Andreas Muller, and is available on the University of Potsdam's TOCSY -
% Toolboxes for Complex Systems: http://tocsy.pik-potsdam.de/petropy.php

function H = petropy(x,n,tau,method,accu)
% The petropy function calculates the permutation entropy of data series.
%
% version 1.1, 27.02.2015:
%   - corrected typo in the description ('same' -> 'equal')
%   - line 106: changed unique-function in newer MATLAB-version
%
% Permutation Entropy
%
% H = PETROPY(X,N,TAU,METHOD,ACCU) computes the permutation entropy H of
% a scalar vector X, using permutation order N, time lags from TAU and
% METHOD to treat equal values. The ACCU parameter describes the accuracy
% of the values in X by the number of decimal places.
%
% x      - data vector (Mx1 or 1xM)
% n      - permutation order
% tau    - time lag scalar OR time lag vector (length = n-1)
% method - method how to treat equal values
%   'noise' - add small noise
%   'equal' - allow same rank for equal values
%   'order' - consider order of appearance (first occurence --> lower rank)
% accu   - maximum number of decimal places in x
%         (only used for method 'noise')
%
% References:
%
% Bandt, C.; Pompe, B. Permutation Entropy: A Natural Complexity
% Measure for  Time Series. Phys. Rev. Lett. 88 (2002) 17, 174102
%
% Riedl, M.; Müller, A.; Wessel, N.: Practical considerations of
% permutation entropy. The European Physical Journal Special Topics
% 222 (2013) 2, 249–262
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example:
%
% H = petropy([6,9,11,12,8,13,5],3,1,'order');
% H =
%       1.5219
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5, accu = 4; end
if nargin < 4, method = 'order'; end

x = x(:);
M = length(x);
equal = false;

if n*log10(n)>15, error('permutation dimension too high'); end
if (length(tau)) > 1 && (length(tau) ~= n-1), error('time lag vector has to have n-1 entries'); end
if ((n-1)*min(tau) >=M) || max(tau) >= M, error('too few data points for desired dimension and lags'); end


switch lower(method)
    case {'noise'}
        %disp('Method: add small noise')
        x = x + rand(M,1)*10^(-accu-1);
    case 'equal'
        %disp('Method: allow equal ranks')
        equal = true;
    case 'order'
        %disp('Method: consider order of occurrence')
    otherwise
        error('unknown method')
end

if length(tau) > 1
    tau = reshape(tau,length(tau),1);
    tau = sort([0;tau]);
    % build n x (M-tau(n))-matrix from shifted values in x
    shift_mat = zeros(n,M-tau(n));
    for ii=1:n
        shift_mat(ii,:) = x(tau(ii)+1:M-tau(n)+tau(ii));
    end
else
    % vectorized
    shift_mat_ind = reshape(0:tau:(n-1)*tau,[],1) * ones(1,M-(n-1)*tau) +...
        ones(n, 1) * reshape(1:(M-(n-1)*tau),1,[]);
    shift_mat = x(shift_mat_ind);
end

if equal
    % allow equal values the same index
    ind_mat = zeros(size(shift_mat));
    for ii=1:size(ind_mat,2)
        [~,~,ind_mat(:,ii)]=unique(shift_mat(:,ii),'first');
    end
else
    % sort matrix along rows to build rank orders, equal values retain
    % order of appearance
    [~, sort_ind_mat] = sort(shift_mat,1);
    ind_mat = zeros(size(sort_ind_mat));
    for ii=1:size(ind_mat,2)
        ind_mat(sort_ind_mat(:,ii),ii) = 1:n;
    end
end
% assign unique number to each pattern (base-n number system)
ind_vec = n.^(0:n-1) * (ind_mat-1);

% find first occurence of unique values in 'ind_vec' and use
% difference to determine length of sequence of the same numbers; e.g.
% sort_ind_vec = [21 21 11 19 11], unique_values = [11 19 21],
% ia = [1 3 4]: 11 occurs on places #1 and #2, 19 on #3 and 21 on #4 and #5
[~,ia,~] = unique(sort(ind_vec), 'first');

% use the following line, if you are using MATLAB in a lower version than
% 8.3:
% permpat_num = diff([ia (length(ind_vec)+1)]);
permpat_num = diff([ia; (length(ind_vec)+1)]);

permpat_num = permpat_num/sum(permpat_num);

% compute permutation entropy
H = -sum(permpat_num .* log2(permpat_num));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code for the 0-1 test is modified from code originally
% written by Paul Matthews and made available here:
% https://www.mathworks.com/matlabcentral/fileexchange/25050-0-1-test-for-chaos
%
% The modification is based on the work of Dawes and Freeland (2008), who
% proposed adding a noise term to M(n). The noise is scaled by the variable
% "sig". This modification improves distinguishability between chaotic vs
% strange non-chaotic systems

function [kmedian,p,q]=z1test(x,sig)

if nargin<2
    sig=1;
end

s=size(x);
if s(2)==1
    x=x';
end
N=length(x); j=[1:N];
t=[1:round(N/10)];
M=zeros(1,round(N/10));
c=pi/5+rand(1,100)*3*pi/5;
for its=1:100
    p=cumsum(x.*cos(j*c(its)));q=cumsum(x.*sin(j*c(its)));
    for n=1:round(N/10)
        M(n)=mean( (p(n+1:N)-p(1:N-n)).^2 + (q(n+1:N)-q(1:N-n)).^2 )- ...
            mean(x)^2*(1-cos(n*c(its)))/(1-cos(c(its)))+sig*(rand-.5);
    end
    kcorr(its)=corr(t',M');
end
kmedian=median(kcorr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following implementation of Schreiber's noise reduction algorithm is
% written by Alexandros Leontitsis, as part of the Matlb Chaotic Systems
% Toolbox

function xr=noiserSchreiber(x,K,L,r,repeat,auto)
%Syntax: xr=noiserSchreiber(x,K,L,r,repeat,auto)
%_______________________________________________
%
% Geometrical noise reduction for a time series, accordiong to
% the extremely simple noise reduction method introduced by Screiber.
%
% xr is the vector/matrix with the cleaned time series.
% x is the time series.
% K is the number of dimensions before the corrected point.
% L is the number of dimensions after the corrected point.
% r is the range of the neighborhood.
% repeat is the number of iterations of the whole algorithm.
% auto automaticaly adjusts the range for every iteration if it is set to
%   'auto'.
%
%
% Note:
%
% K+L+1 should be equal to the embedding dimension.
%
%
% Reference:
%
% Schreiber T (1993): Extremely simple noise-reduction method.
% Physical Review E 47: 2401-2404
%
%
% Alexandros Leontitsis
% Department of Education
% University of Ioannina
% Ioannina
% Greece
% e-mail: leoaleq@yahoo.com
% Homepage: http://www.geocities.com/CapeCanaveral/Lab/1421
%
% 14 Jul 2001

if nargin<1 | isempty(x)==1
    error('You should provide a time series.');
else
    % x must be a vector
    if min(size(x))>1
        error('Invalid time series.');
    end
    x=x(:);
    % n is the time series length
    n=length(x);
end

if nargin<2 | isempty(K)==1
    K=1;
else
    % K must be either a scalar or a vector
    if min(size(K))>1
        error('K must be a scalar or a vector.');
    end
    % K must contain integers
    if sum(abs(round(K)-K))~=0
        error('K must contain integers.');
    end
    % K values must be above 1
    if any(K<1)==1
        error('K values must be above 1.');
    end
end

if nargin<3 | isempty(L)==1
    L=1;
else
    % L must be either a scalar or a vector
    if min(size(L))>1
        error('L must be a scalar or a vector.');
    end
    % L must contain integers
    if sum(abs(round(L)-L))~=0
        error('L must contain integers.');
    end
    % L values must be above 1
    if any(L<1)==1
        error('L values must be above 1.');
    end
end

if nargin<4 | isempty(r)==1
    r=std(x);
else
    % r must be either a scalar or a vector
    if min(size(r))>1
        error('r must be a scalar or a vector.');
    end
    % r values must be above 0
    if any(r<0)==1
        error('r values must be above 0');
    end
end

if nargin<5 | isempty(repeat)==1
    repeat=1;
    repeat1=1;
else
    % repeat must be either a scalar or a vector
    if min(size(repeat))>1
        error('repeat must be a scalar or a vector.');
    end
    % repeat must be above 1
    if repeat<1
        error('repeat must be above 1.')
    end
    % repeat must be integer
    if sum(abs(round(repeat)-repeat))~=0
        error('repeat must contain integers.')
    end
    % The elements of repeat must be in increasing order
    if any(diff(repeat)<=0)==1
        error('The elements of repeat must be in increasing order.')
    end
    repeat1=repeat;
    repeat=1:max(repeat);
end

% Only one of K, L, r, or repeat should be vector
l=[length(K),length(L),length(r),length(repeat)];
if length(find(l>1))>1
    error('Only one of dim, tau, r, p, or repeat should be vector.');
end

m=max(l);
K=ones(1,m).*K;
L=ones(1,m).*L;
r=ones(1,m).*r;
repeat=ones(1,m).*repeat;

for i=1:m
    
    if repeat(i)==1
        % Make the phase-space
        [Y,T]=phasespace(x,K(i)+L(i)+1,1);
    else
        % The new Y
        Yr=[];
        [Y,T]=phasespace(xr(1+K(i)*(repeat(i)-1):n-L(i)*(repeat(i)-1),i-1),...
            K(i)+L(i)+1,1);
    end
    
    for j=1:T
        y=Y(j,:);
        lock=radnearest(y,Y,T,r(i),inf);
        lock(find(lock==j))=[];
        if length(lock)==0
            Yr(j,:)=y;
        else
            Ynearest=Y(lock,:);
            if length(lock)==1
                Yr(j,:)=Ynearest;
            else
                Yr(j,:)=mean(Ynearest);
            end
        end
    end
    
    
    % Calculate the reconstructed time series
    xr(:,i)=[zeros(K(i)*repeat(i),1);Yr(:,K(i)+1);zeros(L(i)*repeat(i),1)];
    
    if nargin==6 & auto=='auto'
        j=1+K(i)*repeat(i):n-L(i)*repeat(i);
        if repeat(i)==1
            r(i+1)=std(x(j)-xr(j),1);
        else
            r(i+1)=std(xr(j,i-1)-xr(j,i),1);
        end
    end
end
end

function lock=radnearest(y,Y,T,r,p)
%Syntax: lock=radnearest(y,Y,T,r,p)
%__________________________________
%
% Locks the nearest neighbors of a reference point that lie within a
% radius in a phase-space.
%
% lock returns the points located.
% y is the reference vector.
% Y is the phase space.
% T is the length of the phase space
% r is the radius.
% p defines the norm.
%
% Alexandros Leontitsis
% Department of Education
% University of Ioannina
% 45110 - Dourouti
% Ioannina
% Greece
%
% University e-mail: me00743@cc.uoi.gr
% Lifetime e-mail: leoaleq@yahoo.com
% Homepage: http://www.geocities.com/CapeCanaveral/Lab/1421
%
% June 15, 2001.
% Initialize j
j=0;
% For every phase-space point
for i=1:T
    % Calculate the distance from the reference point
    dist=norm(y-Y(i,:),p);
    % If it is less than r, count it
    if dist<=r
        j=j+1;
        lock(j)=i;
    end
end
if j==0
    lock=NaN;
end
end

function [Y,T]=phasespace(x,dim,tau)
%Syntax: [Y,T]=phasespace(x,dim,tau)
%___________________________________
%
% The phase space reconstruction of a time series x whith the Method Of Delays
% (MOD), in embedding dimension m and for time dalay tau.
%
% Y is the trajectory matrix in the reconstructed phase space.
% T is the phase space length.
% x is the time series.
% dim is the embedding dimension.
% tau is the time delay.
%
%
% Reference:
% Takens F (1981): Detecting strange attractors in turbulence. Lecture notes in
% Mathematics, 898. Springer.
%
%
% Alexandros Leontitsis
% Department of Education
% University of Ioannina
% 45110 - Dourouti
% Ioannina
% Greece
%
% University e-mail: me00743@cc.uoi.gr
% Lifetime e-mail: leoaleq@yahoo.com
% Homepage: http://www.geocities.com/CapeCanaveral/Lab/1421
%
% 11 Mar 2001.
if nargin<1 | isempty(x)==1
    error('You should provide a time series.');
else
    % x must be a vector
    if min(size(x))>1
        error('Invalid time series.');
    end
    x=x(:);
    % N is the time series length
    N=length(x);
end
if nargin<2 | isempty(dim)==1
    dim=2;
else
    % dim must be scalar
    if sum(size(dim))>2
        error('dim must be scalar.');
    end
    % dim must be an integer
    if dim-round(dim)~=0
        error('dim must be an integer.');
    end
    % dim must be positive
    if dim<=0
        error('dim must be positive.');
    end
end
if nargin<3 | isempty(tau)==1
    tau=1;
else
    % tau must be scalar
    if sum(size(tau))>2
        error('tau must be scalar.');
    end
    % tau must be an integer
    if tau-round(tau)~=0
        error('tau must be an integer.');
    end
    % tau must be positive
    if tau<=0
        error('tau must be positive.');
    end
end
% Total points on phase space
T=N-(dim-1)*tau;
% Initialize the phase space
Y=zeros(T,dim);
% Phase space reconstruction with MOD
for i=1:T
    Y(i,:)=x(i+(0:dim-1)*tau)';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following surrogate data generator is provided by G. Lancaster and
% colleagues. See reference and distrubtion rights below


% Version 1.0 20/06/2018
% *************************************************************************
% ********************* Surrogate data generator **************************
% *************************************************************************
% -----------------------------Copyright-----------------------------------
%
% Authors: Dmytro Iatsenko & Gemma Lancaster
% This software accompanies the article "Surrogate data for hypothesis
% testing in physical systems", G. Lancaster, D. Iatsenko, A. Pidde,
% V. Ticcinelli and A. Stefanovska. Physics Reports, 2018.
%
% Bug reports, suggestions and comments are welcome. Please email them to
% physics-biomed@lancaster.ac.uk
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% surrogate.m is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% See <http://www.gnu.org/licenses/>.
%
%-----------------------------Documentation--------------------------------
% EXAMPLE [surr,params]=surrogate(x, 100, 'FT', 1, 40)
% calculates 100 FT surrogates of the input signal x, which is sampled at
% 40Hz, with preprocessing on.
%
% AVAILABLE SURROGATES:
% 'RP' - Random permutation
% 'FT' - Fourier transform (see also: J. Theiler, S. Eubank, A. Longtin,
% B. Galdrikian, J. Farmer, Testing for nonlinearity in time series: The
% method of surrogate data, Physica D 58 (1�4) 13 (1992) 77�94).
% 'AAFT' - Amplitude adjusted Fourier transform
% 'IAAFT1' - Iterative amplitude adjusted Fourier transform with exact
% distribution (see also: T. Schreiber, A. Schmitz, Improved surrogate data for
% nonlinearity tests, Phys. Rev. Lett. 77 (4) (1996) 635�638).
% 'IAAFT2' - Iterative amplitude adjusted Fourier transform with exact
% spectrum
% 'CPP' - Cyclic phase permutation
% 'PPS' - Pseudo-periodic (see also: M. Small, D. Yu, R.G. Harrison, Surrogate
% test for pseudoperiodic time series data, Phys. Rev. Lett. 87 (18) (2001) 188101.)
% 'TS' - Twin (see also: M. Thiel, M.C. Romano, J. Kurths, M. Rolfs, R. Kliegl,
% Twin surrogates to test for complex synchronisation, Europhys. Lett. 75 (4) (2006) 535).
% 'tshift' - Time shifted
% 'CSS' - Cycle shuffled surrogates. Require that the signal can be
% separated into distinct cycles. May require adjustment of peak finding
% parameters.  (see also: J. Theiler, On the evidence for low-dimensional chaos in an epileptic
% electroencephalogram, Phys. Lett. A 196 (1) (1995) 335�341).
%
%
% INPUT:
% sig - time series for which to calculate surrogate(s)
% Input sig should be the original time series, or phase for CPP.
% For surrogates requiring embedding, the delay tau and dimension D are
% calculated automatically using false nearest neighbours and first 0
% autocorrelation, respectively.
% N - number of surrogates to calculate
% method - surrogate type, choose from one of the strings above
% pp - preprocessing on (1) or off (0) (match beginning and end and first
% derivatives)
%
% varargin options
% if method = 'FT', random phases can be input and output, to preserve for
% multivariate surrogates for example
% if method = 'PPS' or 'TS', embedding dimension can be entered beforehand instead of
% being estimated, if it is to be kept the same for all surrogates
% if method = 'CSS', minimum peak height and minimum peak distance can be
% entered to ensure correct peak detection for separation of cycles.
%
% OUTPUT:
% surr - surrogate data
% params - parameters from the generation, including truncation locations
% if preprocessed and runtime.

function [surr,params]=surrogate(sig, N, method, pp, fs, varargin)

origsig=sig;
params.origsig=origsig;
params.method=method;
params.numsurr=N;
params.fs=fs;
z=clock;

%%%%%%% Preprocessing %%%%%%%
%%
if pp==1
    [sig,time,ks,ke]=preprocessing(sig,fs);
else
    time=linspace(0, length(sig)/fs,length(sig));
end
L=length(sig);
L2=ceil(L/2);
if pp==1
    params.preprocessing='on';
    params.cutsig=sig;
    params.sigstart=ks;
    params.sigend=ke;
else
    params.preprocessing='off';
end
params.time=time;
%%
%%%%%%% Random permutation (RP) surrogates %%%%%%%
if strcmp(method,'RP')
    surr=zeros(N,length(sig));
    for k=1:N
        surr(k,:)=sig(randperm(L));
    end
    
    
    %%%%%%% Fourier transform (FT) surrogate %%%%%%%
elseif strcmp(method,'FT')
    
    a=0; b=2*pi;
    if nargin>5
        eta=varargin{1};
    else
        eta=(b-a).*rand(N,L2-1)+a; % Random phases
    end
    ftsig=fft(sig); % Fourier transform of signal
    ftrp=zeros(N,length(ftsig));
    
    ftrp(:,1)=ftsig(1);
    F=ftsig(2:L2);
    F=F(ones(1,N),:);
    ftrp(:,2:L2)=F.*(exp(1i*eta));
    ftrp(:,2+L-L2:L)=conj(fliplr(ftrp(:,2:L2)));
    
    surr=ifft(ftrp,[],2);
    
    params.rphases=eta;
    
    
    
    %%%%%%% Amplitude adjusted Fourier transform surrogate %%%%%%%
elseif strcmp(method,'AAFT')
    
    a=0; b=2*pi;
    eta=(b-a).*rand(N,L2-1)+a; % Random phases
    [val,ind]=sort(sig);
    rankind(ind)=1:L;    % Rank the locations
    
    gn=sort(randn(N,length(sig)),2); % Create Gaussian noise signal and sort
    for j=1:N
        gn(j,:)=gn(j,rankind); % Reorder noise signal to match ranks in original signal
    end
    
    ftgn=fft(gn,[],2);
    F=ftgn(:,2:L2);
    
    surr=zeros(N,length(sig));
    surr(:,1)=gn(:,1);
    surr(:,2:L2)=F.*exp(1i*eta);
    surr(:,2+L-L2:L)=conj(fliplr(surr(:,2:L2)));
    surr=(ifft(surr,[],2));
    
    [~,ind2]=sort(surr,2); % Sort surrogate
    rrank=zeros(1,L);
    for k=1:N
        rrank(ind2(k,:))=1:L;
        surr(k,:)=val(rrank);
    end
    
    
    
    
    %%%%%%% Iterated amplitude adjusted Fourier transform (IAAFT-1) with exact distribution %%%%%%%
elseif strcmp(method,'IAAFT1')
    maxit=1000;
    [val,ind]=sort(sig);  % Sorted list of values
    rankind(ind)=1:L; % Rank the values
    
    ftsig=fft(sig);
    F=ftsig(ones(1,N),:);
    surr=zeros(N,L);
    
    for j=1:N
        surr(j,:)=sig(randperm(L)); % Random shuffle of the data
    end
    
    it=1;
    irank=rankind;
    irank=irank(ones(1,N),:);
    irank2=zeros(1,L);
    oldrank=zeros(N,L);
    iind=zeros(N,L);
    iterf=zeros(N,L);
    
    while max(max(abs(oldrank-irank),[],2))~=0 && it<maxit
        go=max(abs(oldrank-irank),[],2);
        [~,inc]=find(go'~=0);
        
        oldrank=irank;
        iterf(inc,:)=real(ifft(abs(F(inc,:)).*exp(1i*angle(fft(surr(inc,:),[],2))),[],2));
        
        [~,iind(inc,:)]=sort(iterf(inc,:),2);
        for k=inc
            irank2(iind(k,:))=1:L;
            irank(k,:)=irank2;
            surr(k,:)=val(irank2);
        end
        
        it=it+1;
    end
    
    
    %%%%%%% Iterated amplitude adjusted Fourier transform (IAAFT-2) with exact spectrum %%%%%%%
elseif strcmp(method,'IAAFT2')
    maxit=1000;
    [val,ind]=sort(sig);  % Sorted list of values
    rankind(ind)=1:L; % Rank the values
    
    ftsig=fft(sig);
    F=ftsig(ones(1,N),:);
    surr=zeros(N,L);
    
    for j=1:N
        surr(j,:)=sig(randperm(L)); % Random shuffle of the data
    end
    
    it=1;
    irank=rankind;
    irank=irank(ones(1,N),:);
    irank2=zeros(1,L);
    oldrank=zeros(N,L);
    iind=zeros(N,L);
    iterf=zeros(N,L);
    
    while max(max(abs(oldrank-irank),[],2))~=0 && it<maxit
        go=max(abs(oldrank-irank),[],2);
        [~,inc]=find(go'~=0);
        
        oldrank=irank;
        iterf(inc,:)=real(ifft(abs(F(inc,:)).*exp(1i*angle(fft(surr(inc,:),[],2))),[],2));
        
        [~,iind(inc,:)]=sort(iterf(inc,:),2);
        for k=inc
            irank2(iind(k,:))=1:L;
            irank(k,:)=irank2;
            surr(k,:)=val(irank2);
        end
        
        it=it+1;
        
    end
    surr=iterf;
    
    
    
    %%%%%%% Cyclic phase permutation (CPP) surrogates %%%%%%%
elseif strcmp(method,'CPP')
    
    phi=wrapTo2Pi(sig);
    
    pdiff=phi(2:end)-phi(1:end-1);
    locs=find(pdiff<-pi);
    if length(locs)==1
        locs=[1 locs];
    end 
    parts=cell(length(locs)-1);
    for j=1:length(locs)-1
        tsig=phi(locs(j)+1:locs(j+1));
        parts{j}=tsig;
    end
    
    st=phi(1:locs(1));
    en=phi(locs(j+1)+1:end);
    surr=zeros(N,L);
    for k=1:N
        surr(k,:)=unwrap(horzcat(st,parts{randperm(j)},en));
        
    end
    
    
    %%%%%%% Pseudo-periodic surrogates (PPS) %%%%%%%
elseif strcmp(method,'PPS')
    
    % Embedding of original signal
    if nargin>5
        m=varargin{1};
        [sig,tau]=embedsig(sig,'DimAlg',m);
        
    else
        [sig,tau]=embedsig(sig,'DimAlg','fnn');
        m=size(sig);
        m=m(1);
    end
    
    L=length(sig);
    L2=ceil(L/2);
    time=linspace(0,length(sig)/fs,length(sig));
    params.embed_delay=tau;
    params.embed_dim=m;
    params.embed_sig=sig;
    
    % % Find the index of the first nearest neighbour from the first half of the
    % % embedded signal to its last value to avoid cycling near last value
    
    for k=1:L
        matr=max(abs(sig(:,:)-sig(:,k)*ones(1,L)));
        [ssig(k),mind(k)]=min(matr(matr>0));
    end
    
    [~,pl]=min(matr(1:round(L/2))); rho=0.7*mean(ssig); clear mind ssig;
    parfor x=1:N
        kn=randi(L,1); % Choose random starting point
        
        for j=1:L % Length of surrogate is the same as the embedded time series
            if kn==L
                kn=pl;
            end
            kn=kn+1; % Move forward from previous kn
            surr(x,j)=sig(1,kn); % Set surrogate to current value for kn (choose first component, can be any)
            sigdist=max(abs(sig(:,:)-(sig(:,kn)+randn*rho)*ones(1,L))); % Find the maximum
            % distance between each point in the original signal and the current
            % values with noise added
            [~,kn]=min(sigdist); % Find nearest point
            
        end
    end
    
    
    
    % %%%%%%% Twin surrogates %%%%%%%
elseif strcmp(method,'TS') %
    
    % Embedding of original signal
    if nargin>5
        m=varargin{1};
        [sig,tau]=embedsig(sig,'DimAlg',m);
        
    else
        [sig,tau]=embedsig(sig,'DimAlg','fnn');
        m=size(sig);
        m=m(1);
    end
    L=length(sig);
    L2=ceil(L/2);
    time=linspace(0,length(sig)/fs,length(sig));
    params.embed_delay=tau;
    params.embed_dim=m;
    params.embed_sig=sig;
    
    dL=L;
    alpha=0.1;
    
    Rij=zeros(L,L);
    for k=2:L
        Rij(k,1:k-1)=max(abs(sig(:,1:k-1)-sig(:,k)*ones(1,k-1)));
    end
    Rij=Rij+Rij';
    [~,pl]=min(Rij(1:round(L/2),L));
    Sij=sort(Rij(:)); delta=Sij(round(alpha*L^2)); clear Sij;
    Rij(Rij<delta)=-1; Rij(Rij>delta)=0; Rij=abs(Rij);
    
    ind=cell(L,1); eln=zeros(L,1); twind=1:L;
    remp=1; % remaining points
    while ~isempty(remp)
        twn=remp(1);
        ind{twn}=remp(max(abs(Rij(:,remp)-Rij(:,twn)*ones(1,numel(remp))))==0);
        ind(ind{twn})=ind(twn);
        eln(ind{twn})=length(ind{twn});
        twind(ind{twn})=0;
        remp=twind(twind>0);
    end
    clear Rij twind;
    
    for sn=1:N
        kn=randi(L,1)-1;
        for j=1:dL
            kn=kn+1;
            surr(sn,j)=sig(1,kn);
            kn=ind{kn}(randi(eln(kn),1));
            if kn==L
                kn=pl;
            end
        end
    end
    
    
    
    %%%%%%% Time-shifted surrogates %%%%%%%
elseif strcmp(method,'tshift')
    %nums=randperm(L);
    for sn=1:N
        startp=randi(L-1,1);%nums(sn);%
        surr(sn,:)=horzcat(sig(1+startp:L),sig(1:startp));
    end
    %params.tshifts=nums(1:N);
    
    
    
    
    %%%%%%% Cycle shuffled surrogates
elseif strcmp(method,'CSS')
    
    if nargin>5
        MPH=varargin{1}; % Minimum heak height
        MPD=varargin{2}; % Minimum peak distance
    else
        MPH=0;
        MPD=fs;
    end
    
    [~,I]=findpeaks(sig,'MinPeakHeight',MPH,'MinPeakDistance',MPD);
    
    st=sig(1:I(1)-1);
    en=sig(I(end):end);
    
    for j=1:length(I)-1
        parts{j}=sig(I(j):I(j+1)-1);
    end
    
    for k=1:N
        surr(k,:)=unwrap(horzcat(st,parts{randperm(j)},en));
        
    end
    
end

params.runtime=etime(clock,z);
params.type=method;


end




function [cutsig,t2,kstart,kend]=preprocessing(sig,fs)
sig=sig-mean(sig);
t=linspace(0,length(sig)/fs,length(sig));
L=length(sig);
p=10; % Find pair of points which minimizes mismatch between p consecutive
%points and the beginning and the end of the signal

K1=round(L/100); % Proportion of signal to consider at the beginning
k1=sig(1:K1);
K2=round(L/10);  % Proportion of signal to consider at the end
k2=sig(end-K2:end);

% Truncate to match start and end points and first derivatives
if length(k1)<=p
    p=length(k1)-1;
else
end
d=zeros(length(k1)-p,length(k2)-p);

for j=1:length(k1)-p
    for k=1:length(k2)-p
        d(j,k)=sum(abs(k1(j:j+p)-k2(k:k+p)));
    end
end

[v,I]=min(abs(d),[],2);
[~,I2]=min(v); % Minimum mismatch

kstart=I2;
kend=I(I2)+length(sig(1:end-K2))-1;
cutsig=sig(kstart:kend); % New truncated time series
t2=t(kstart:kend); % Corresponding time

end

function vec = minmaxsig(a)
inds1=islocalmax(a);
inds2=islocalmin(a);
inds=logical(inds1+inds2);
vec=a(inds);
end


% The following is code for Breitung's variance ratio test. Code and lookup
% tables were adapted from the R code written by Matthew Clegg as part
% of the Engle-Granger Cointegration Models (EGCM) toolbox: 
%https://github.com/cran/egcm/tree/master/man

function [H,rho,pval] = bvr(y, dt,alpha)
% Calculates the variance ratio statistic rho described in equation (5) of
%   Breitung, Jorg (2001).  Nonparametric tests for unit roots and cointegration,
%   Journal of Econometrics, 108, 343-363.
if nargin<3
    alpha=0.05;
end
if dt==1
    y = detrend(y); 
    qtab = [NaN NaN 50 100 250 500 750 1000 1250;0.1 0.001...
        0.001408808 0.00133 0.0014 0.00138 0.00137 0.00137 0.00135;...
        1 0.01 0.002299207 0.0022 0.0022 0.0022 0.00221 0.00221 0.0022;...
        2.5 0.025 0.00283679 0.00282 0.00278 0.00278 0.00278 0.0028 0.00279;...
        5 0.05 0.003499243 0.00346 0.00341 0.00341 0.00344 0.00346 0.00343;...
        10 0.1 0.004432426 0.00441 0.00443 0.00438 0.00441 0.00439 0.00437;...
        20 0.2 0.006016803 0.00598 0.00601 0.00591 0.00593 0.00596 0.00595;...
        50 0.5 0.010274684 0.0102 0.0102 0.0101 0.0101 0.0102 0.0101;...
        80 0.8 0.015998227 0.0159 0.0159 0.0159 0.0158 0.016 0.0158;...
        90 0.9 0.018588419 0.0185 0.0185 0.0186 0.0185 0.0185 0.0184;...
        95 0.95 0.020253402 0.0202 0.0201 0.0202 0.0201 0.0202 0.0201;...
        97.5 0.975 0.021314538 0.0213 0.0212 0.0213 0.0212 0.0213 0.0212;...
        99 0.99 0.022280874 0.0222 0.0222 0.0222 0.0222 0.0222 0.0221;...
        99.9 0.999 0.0236688 0.0234 0.0236 0.0235 0.0235 0.0235 0.0234];  
elseif dt==0
    y = y - mean(y);
    qtab=[NaN NaN 50 100 250 500 750 1000 1250;0.1 0.001 0.003129817 0.00303...
        0.0029 0.00306 0.00301 0.00298 0.00306;1 0.01 0.005804221 0.00559...
        0.00553 0.00545 0.00549 0.00551 0.00538;2.5 0.025 0.007848024 0.00772...
        0.00763 0.00755 0.00767 0.00767 0.00747;5 0.05 0.010367278 0.0104...
        0.0101 0.0101 0.0102 0.0103 0.0101;10 0.1 0.014482161 0.0146 0.0143...
        0.0142 0.0144 0.0144 0.0143;20 0.2 0.021578693 0.0217 0.0215 0.0213...
        0.0215 0.0214 0.0216;50 0.5 0.051426484 0.0513 0.0515 0.0512 0.0509...
        0.0512 0.0516;80 0.8 0.079608597 0.0798 0.0794 0.0794 0.0795 0.0794...
        0.0797;90 0.9 0.08743295 0.0876 0.0874 0.0874 0.0876 0.0875 0.0878;...
        95 0.95 0.091724245 0.0919 0.0916 0.0917 0.0919 0.0917 0.0918;97.5...
        0.975 0.094158785 0.0943 0.0942 0.0942 0.0943 0.0941 0.0942;99 0.99...
        0.096180563 0.0963 0.0963 0.096 0.0962 0.0962 0.0963;99.9 0.999...
        0.098546067 0.0986 0.0987 0.0985 0.0987 0.0985 0.0986];
end

ys = cumsum(y);
n = length(y);
rho_num = (1 / n^2)*sum(ys.^2);
rho_den = (1 / n) * sum(y.^2);
% Note factor of 1/n has been added that was not in the original paper.
rho = (rho_num / rho_den) / n;
pval=quantile_table_interpolate(qtab,length(y),rho);
H= pval<=0.05;
    
end


function pval = quantile_table_interpolate(qtab, sample_size, stat)

n = size(qtab,1);
i=0;
for j = 2:size(qtab,2)-1
    if qtab(1,j)<=sample_size && qtab(1,j+1)>sample_size
        i=j;
    end
    if i==0
        i=size(qtab,2);
    end
end
if stat<min(qtab(2:n,i))
    y1=min(qtab(2:n,i));
elseif stat>max(qtab(2:n,i))
    y1=max(qtab(2:n,i));
else
    y1 = interp1(qtab(2:n,i),qtab(2:n,2),stat);
end
if i < size(qtab,2)
    if stat<min(qtab(2:n,i+1))
        y2=min(qtab(2:n,i+1));
    elseif stat>max(qtab(2:n,i+1))
        y2=max(qtab(2:n,i+1));
    else
        y2 = interp1(qtab(2:n,i),qtab(2:n,2),stat);
    end
    n1 = qtab(1,i);
    n2 = qtab(1,i+1);
    pval = y1*(n2 - sample_size) / (n2 - n1) + y2 * (sample_size - n1)/(n2 - n1);
else
    pval=y1;
end
end

function gauss_dat = rank_inverse_gaussian_transform(y)

% rank-based inverse normal transformation
if size(y,1)>size(y,2)
    y=y';
end
rank = tiedrank(y);
p = rank / (length(rank)+1);
gauss_dat = norminv( p, 0, 1 );
end

