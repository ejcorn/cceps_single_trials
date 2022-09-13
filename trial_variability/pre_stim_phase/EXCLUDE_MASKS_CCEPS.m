function [excl_elecs,excl_tps] = EXCLUDE_MASKS_CCEPS(X,good_cceps,trial_len)

% INPUTS:
% X: TxN processed CCEP trial by trial data matrix with some missingness
% due to high amplitude signals most often in or near stim electrodes
% good_cceps: NxN matrix with 1's indicating which CCEPs are suprathreshold
% trial_len: number of time points per trial
%
% OUTPUTS:
% excl_chans: mask with 1's for channels to keep
% excl_chans: mask with 1's for time points to keep


% initialize exclusion masks
excl_elecs = true(1,size(X,2));
excl_tps = true(size(X,1),1);        

% exclude stim electrodes
good_stim_elecs = find(any(good_cceps,1));
excl_elecs(good_stim_elecs) = 0;

% exclude electrodes that are blank in a bipolar montage or missing
% from every single trial
chan_nan_thresh = 0.5; %1/length(X); % exclude channels if >50% nan
excl_elecs(mean(isnan(X),1) >= chan_nan_thresh) = 0;        

% create index for trials
n_trials = size(X,1)/trial_len;
trial_idx = repelem(1:n_trials,trial_len)';        

% exclude bad trials
tp_nan_thresh = 0; %remove any trials with any nans in any electrode
% find trials (trial_idx == x) with any nans in it
bad_trials = cell2mat(cellfun(@(x) mean(isnan(X(trial_idx==x,excl_elecs)),[1 2]) > tp_nan_thresh,num2cell(1:n_trials),'UniformOutput',false))';
% exclude the time points belonging to those trials
excl_tps(ismember(trial_idx,find(bad_trials))) = 0;
