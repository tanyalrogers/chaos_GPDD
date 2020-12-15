function G = HVG_alg(ts,left,right,G,t,sh)
%   Fast Horizontal visibility algorithm  based on:
%   Lan, X., Mo, H., Chen, S., Liu, Q., & Deng, Y. (2015).
%   Fast transformation from time series to visibility graphs.
%   Chaos, 25(8), 083105.
%
%   For the Natural visibility algorithm, see Giovanni Iacobello (2019).
%   Fast natural visibility graph (NVG) for MATLAB, MATLAB Central File Exchange
%   ===============================================================
%   Code by: Giovanni Iacobello, Politecnico di Torino (Italy)
%   giovanni.iacobello@polito.it
%   ===============================================================
%   ts=time series
%   left= first data index
%   right= last data index
%   G=adjacency list (cell array)
%   t=time samples
%   sh=series vertical shifting (real value)

if left<right
    [~,k]=max(ts(left:right));
    k=k+left-1;
    
    beta=-Inf;
    neigs_k=zeros(right-left+1,1);
    cn=1;
    for ii=k-1:-1:left
        alfa=ts(ii);
        if alfa>beta+sh
           neigs_k(cn)=ii;
           beta=alfa;
           cn=cn+1;
        end
    end
    beta=-Inf;
    for ii=k+1:right
        alfa=ts(ii);
        if alfa>beta+sh
          neigs_k(cn)=ii;
           beta=alfa;
           cn=cn+1;
        end
    end
    G{k,1}=neigs_k(1:cn-1); 
    G=HVG_alg(ts,left,k-1,G,t,sh);
    G=HVG_alg(ts,k+1,right,G,t,sh);
end

return;
end   