function [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_SEQ_SPEARMAN(out,data)

% INPUTS:
% data: NxN cell, each element contains Tx1 vector containing waveform or
% other metric for each of T trials (in order)
% A: CCEPs matrix for additional thresholding of N1 and N2 to discard bad
% average waveforms
% 
% OUTPUTS:
% spear_trials: NxN matrix containing spearman between trial order and
% trial waveform/metric
% spear_trials_p: NxN matrix containing p value for above spearman correlation
% spear_trials_p_adj: spear_p, FDR corrected

good_cceps = ~isempty_c(data) & out.network(1).A>0 & out.network(2).A>0;

% exclude cell elements that have too few non-nan data points - specifying
% must be at least 10 trials

non_nan_count = repmat({false},length(data));
non_nan_count(good_cceps) = cellfun(@(x) sum(~isnan(x))>10,data(good_cceps),'UniformOutput',false);
non_nan_count = cell2mat(non_nan_count);

% use that as new good_cceps mask since it is a subset
good_cceps = good_cceps & non_nan_count;

spear = @(x) corr(x,[1:size(x,1)]','rows','pairwise','type','spearman');
spear_trials = repmat({nan},length(data));
spear_trials(good_cceps) = cellfun(spear,data(good_cceps),'UniformOutput',false);
spear_trials = cell2mat(spear_trials);

spear_p = @(x) CORR_SPEAR_P(x,[1:size(x,1)]');
spear_trials_p = repmat({nan},length(data));
spear_trials_p(good_cceps) = cellfun(spear_p,data(good_cceps),'UniformOutput',false);
spear_trials_p = cell2mat(spear_trials_p);
spear_trials_p_adj = nan(size(spear_trials_p));
[~,~,~,spear_trials_p_adj(good_cceps)] = fdr_bh(spear_trials_p(good_cceps));

%{
for j = find(good_cceps)'
    
    spear_p(data{j})
end
%}