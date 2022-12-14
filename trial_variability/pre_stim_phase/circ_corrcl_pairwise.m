function [rho pval] = circ_corrcl_pairwise(alpha, x)
%
% [rho pval ts] = circ_corrcc(alpha, x)
%   Correlation coefficient between one circular and one linear random
%   variable.
%
%   Input:
%     alpha   sample of angles in radians
%     x       sample of linear random variable
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value
%
% References:
%     Biostatistical Analysis, J. H. Zar, p. 651
%
% PHB 6/7/2008
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

% MODIFIED 5/28/22 to ensure only rows that are not nan for both x and
% alpha are included - otherwise using the 'pairwise' rows argument leads
% to inaccurate computation relying on different observations for rxs/rxc
% vs. rcs (lines 51-53 below)


if size(alpha,2) > size(alpha,1)
	alpha = alpha';
end

if size(x,2) > size(x,1)
	x = x';
end

if length(alpha)~=length(x)
  error('Input dimensions do not match.')
end

n = length(alpha);

% added by EJC - 
mask = ~isnan(alpha) & ~isnan(x);
x = x(mask);
alpha = alpha(mask);

% compute correlation coefficent for sin and cos independently
rxs = corr(x,sin(alpha));
rxc = corr(x,cos(alpha));
rcs = corr(sin(alpha),cos(alpha));

% compute angular-linear correlation (equ. 27.47)
rho = sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1-rcs^2));

% compute pvalue
pval = 1 - chi2cdf(n*rho^2,2);

