function fail_rate = COMPUTE_FAILURE_RATE(out,data)

% INPUTS:
% out: CCEPs out struct for corresponding patient
% data: NxN cell, each element contains Tx1 vector containing waveform or
% other metric for each of T trials (in order)
% 
% OUTPUTS:
% fail_rate: NxN matrix containing failure rate (i.e. percent of
% good,non-artifactual trials for each CCEP that have undetectable peaks 
% i.e. exactly == 0

good_cceps = ~isempty_c(data) & out.network(1).A>0 & out.network(2).A>0;

% exclude cell elements that have too few non-nan data points - specifying
% must be at least 10 trials

non_nan_count = repmat({false},length(data));
non_nan_count(good_cceps) = cellfun(@(x) sum(~isnan(x))>10,data(good_cceps),'UniformOutput',false);
non_nan_count = cell2mat(non_nan_count);

% use that as new good_cceps mask since it is a subset
good_cceps = good_cceps & non_nan_count;

fail_rate = repmat({nan},length(data));
fail_rate(good_cceps) = cellfun(@(x) nanmean(x==0),data(good_cceps),'UniformOutput',false);
fail_rate = cell2mat(fail_rate);
