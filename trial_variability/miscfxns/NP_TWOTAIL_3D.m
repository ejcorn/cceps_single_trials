function p = NP_TWOTAIL_3D(dist,test_val)
% INPUTS:
% test_val: individual value being compared to distribution
% dist: KxRxN matrix,distribution of values under some null model, or otherwise
% K and R are just dimensions of matrix of variables, N is number of
% observations for each variable
%
% compute 2-tailed p-value for test value occurring in distribution
%
% OUTPUTS:
% p: KxRx1 vector of two-tailed p-values

p = 2*min(cat(3,mean(dist <= test_val,3),mean(dist >= test_val,3)),[],3);