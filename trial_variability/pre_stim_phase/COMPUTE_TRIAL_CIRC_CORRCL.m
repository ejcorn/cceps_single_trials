function [clcorr_trials,clcorr_trials_p,clcorr_trials_p_adj] = COMPUTE_TRIAL_CIRC_CORRCL(out,data_circ,data_lin)

% INPUTS:
% out: processed data structure for patient
% data_circ: NxN cell, each element contains Tx1 vector containing phase data or
% other CIRCULAR metric for each of T trials (in order)
% data_lin: NxN cell, each element contains Tx1 vector containing waveform or
% other LINEAR metric for each of T trials (in order)
% 
% OUTPUTS:
% clcorr_trials: NxN matrix containing circular-linear correlation between trial order and
% trial waveform/metric
% clcorr_trials_p: NxN matrix containing p value for above clcorrman correlation
% clcorr_trials_p_adj: clcorr_p, FDR corrected

good_cceps = ~isempty_c(data_circ) & ~isempty_c(data_lin) & out.network(1).A>0 & out.network(2).A>0;

if ~isequal(size(data_circ),size(data_lin))
    error('Data cells are of unequal sizes');
end

% exclude cell elements that have too few non-nan data points (i.e. not
% enough to get a p value from t test so 3
non_nan_count = repmat({false},length(data_circ));
non_nan_count(good_cceps) = cellfun(@(x,y) sum(~isnan(x) & ~isnan(y))>10,data_circ(good_cceps),data_lin(good_cceps),'UniformOutput',false);
non_nan_count = cell2mat(non_nan_count);

% use that as new good_cceps mask since it is a subset
good_cceps = good_cceps & non_nan_count;

clcorr = @(alpha,x) circ_corrcl_pairwise(alpha,x);
clcorr_trials = repmat({nan},length(data_circ));
clcorr_trials(good_cceps) = cellfun(clcorr,data_circ(good_cceps),data_lin(good_cceps),'UniformOutput',false);
clcorr_trials = cell2mat(clcorr_trials);

clcorr_p = @(alpha,x) CIRC_CORRCL_P(alpha,x);
clcorr_trials_p = repmat({nan},length(data_circ));
clcorr_trials_p(good_cceps) = cellfun(clcorr_p,data_circ(good_cceps),data_lin(good_cceps),'UniformOutput',false);
clcorr_trials_p = cell2mat(clcorr_trials_p);
clcorr_trials_p_adj = nan(size(clcorr_trials_p));
[~,~,~,clcorr_trials_p_adj(good_cceps)] = fdr_bh(clcorr_trials_p(good_cceps));
    