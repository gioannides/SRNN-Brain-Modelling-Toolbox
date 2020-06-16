%-----------------------------------------------------------------------
% FUNCTION: aks_diff.m
% PURPOSE:  apply differencing to a data matrix
% 
% INPUTS:   M:      nvar x nobs data matrix
%               
% OUTPUT:   M2:     differenced data matrix (nvar x nobs-1)   
%     
%           Written by Anil Seth, December 2005
%           Ref: Seth, A.K. (2005) Network: Comp. Neural. Sys. 16(1):35-55
%-----------------------------------------------------------------------
function [M2] = aks_diff(M)

[nvar nobs] = size(M);
if(nobs < nvar) 
    error('error in aks_diff: nobs < nvar, check data matrix');
end
M2 = diff(M');
M2 = M2';
        





