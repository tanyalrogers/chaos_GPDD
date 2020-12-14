function [VG] = fast_HVG(TS,TT,shiftvalue)
% Implementation of horizontal visibility graph based on:
% Luque, B., Lacasa, L., Ballesteros, F., & Luque, J. (2009). 
% Horizontal visibility graphs: Exact results for random time series. 
% Physical Review E, 80(4), 046103.
% ===============================================================
% Code by: Giovanni Iacobello, Politecnico di Torino (Italy)
% giovanni.iacobello@polito.it
% ===============================================================
% TS: time series (must be a vector)
% TT: time vector (same length od TS)
% shiftvalue: scalar real value. Default: shiftvalue=0; 
%             Nodes i and j are linked if
%             (TS(k)+shiftvalue)<min[TS(i),TS(j)], for any i<k<j
% VG: sparse adjacency matrix (symmetrical)
%
%%  PRE-PROCESSING CONTROLS
    VG=[];
    if isvector(TS)==0 || isvector(TT)==0
        disp('Error size: series and times must be vectors')
        return;
    end
    if length(TS)~=length(TT)
        disp('Error size: series and times must have the same length')
        return;
    end
    if iscolumn(TS)==0
        TS=TS';
    end
    if iscolumn(TT)==0
        TT=TT';
    end
    if exist('shiftvalue','var')==0
        shiftvalue=0; %no shifting -> classical HVG
    elseif isscalar(shiftvalue)==0
        disp('Error shiftvalue: insert a scalar value')
        return;
    end    
    
%% RUNNING
    N=size(TS,1);
    B=cell(N,1);
    B=HVG_alg(TS,1,N,B,TT,shiftvalue);
    
%% CONVERTING INTO SPARSE MATRIX

    Adj_tmp=cell(N,1);
    for ii=1:N
        Adj_tmp{ii}=ones(length(B{ii,1}),1)*ii;
    end
    target=vertcat(Adj_tmp{:});
    clear A_tmp
    source=vertcat(B{:,1});
    weight_input=ones(size(target));
    
    if ~isa(source,'double')
        source=double(source); %"sparse" function requires doubles
    end
    VG=sparse(source,target,weight_input,N,N);
    VG=VG+VG'; %makes adjcency matrix simmetrical

end

