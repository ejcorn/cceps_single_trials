function [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_TRIAL_SPEARMAN(out,data1,data2)

% INPUTS:
% out: processed data structure for patient
% data1, data2: NxN cell, each element contains Tx1 vector containing waveform or
% other metric for each of T trials (in order)
% 
% OUTPUTS:
% spear_trials: NxN matrix containing spearman between trial order and
% trial waveform/metric
% spear_trials_p: NxN matrix containing p value for above spearman correlation
% spear_trials_p_adj: spear_p, FDR corrected

good_cceps = ~isempty_c(data1) & ~isempty_c(data2) & out.network(1).A>0 & out.network(2).A>0;

if ~isequal(size(data1),size(data2))
    error('Data cells are of unequal sizes');
end

% exclude cell elements that have too few non-nan data points (i.e. not
% enough to get a p value from t test so 3
non_nan_count = repmat({false},length(data1));
non_nan_count(good_cceps) = cellfun(@(x,y) sum(~isnan(x) & ~isnan(y))>10,data1(good_cceps),data2(good_cceps),'UniformOutput',false);
non_nan_count = cell2mat(non_nan_count);

% use that as new good_cceps mask since it is a subset
good_cceps = good_cceps & non_nan_count;

%spear = @(x,y) corr(x(1:nt),y(1:nt),'rows','pairwise','type','spearman');
spear = @(x,y) corr(x,y,'rows','pairwise','type','spearman');
spear_trials = repmat({nan},length(data1));
spear_trials(good_cceps) = cellfun(spear,data1(good_cceps),data2(good_cceps),'UniformOutput',false);
spear_trials = cell2mat(spear_trials);

%spear_p = @(x,y) CORR_SPEAR_P(x(1:nt),y(1:nt));
spear_p = @(x,y) CORR_SPEAR_P(x,y);
spear_trials_p = repmat({nan},length(data1));
spear_trials_p(good_cceps) = cellfun(spear_p,data1(good_cceps),data2(good_cceps),'UniformOutput',false);
spear_trials_p = cell2mat(spear_trials_p);
spear_trials_p_adj = nan(size(spear_trials_p));
[~,~,~,spear_trials_p_adj(good_cceps)] = fdr_bh(spear_trials_p(good_cceps));
    