function [coefvar_trials,std_trials] = COMPUTE_CROSS_TRIAL_COEFVAR(out,data)

% INPUTS:
% out: CCEPs out struct for corresponding patient
% data: NxN cell, each element contains Tx1 vector containing waveform or
% other metric for each of T trials (in order)
% 
% OUTPUTS:
% std_trials: NxN matrix containing standard deviation of metric across trials ignoring NaNs 

good_cceps = ~isempty_c(data) & out.network(1).A>0 & out.network(2).A>0;

% exclude cell elements that have too few non-nan data points - specifying
% must be at least 10 trials

non_nan_count = repmat({false},length(data));
non_nan_count(good_cceps) = cellfun(@(x) sum(~isnan(x))>10,data(good_cceps),'UniformOutput',false);
non_nan_count = cell2mat(non_nan_count);

% use that as new good_cceps mask since it is a subset
good_cceps = good_cceps & non_nan_count;

coefvar_trials = repmat({nan},length(data));
coefvar_trials(good_cceps) = cellfun(@(x) nanstd(x)/nanmean(x),data(good_cceps),'UniformOutput',false);
coefvar_trials = cell2mat(coefvar_trials);

std_trials = repmat({nan},length(data));
std_trials(good_cceps) = cellfun(@(x) nanstd(x),data(good_cceps),'UniformOutput',false);
std_trials = cell2mat(std_trials);