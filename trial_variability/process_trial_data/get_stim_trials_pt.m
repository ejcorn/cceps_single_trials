function [trials,artifacts,info,trials_raw] = get_stim_trials_pt(out,opts)

% INPUTS:
% out: out struct
% opt: options struct
%
% OUTPUTS:
% trials: NxN cell for N electrodes with time series for each trial
% recorded from each stim-recording pair
%
% if > opt.mt_thresh trials are missing, that stim-rec pair will be
% empty in trials and ignored
%
% artifacts: indices of artifactual trials as determined by artifact
% rejection function
% info: information about each stim-recording pair (index, names, waveform
% strngth), sorted by specified waveform strength

%% Parameters

which = 1; % 1-> N1, 2-> N2
mt_thresh = 0.33; % threshold for number of missing trials

if exist('opts','var')    
    if isfield(opts,'which'); which = opts.which; end
    if isfield(opts,'mt_thresh'); mt_thresh = opts.mt_thresh; end
end

%% Get network
A = out.network(which).A;
nchs = length(A);
stim_idx = repmat(1:size(A,2),size(A,1),1);
response_idx = repmat((1:size(A,1))',1,size(A,2));

%% Find which stim electrodes produced any ccep responses

good_stim_elecs = find(mean(isnan(A))<1);
%good_stim_elecs = good_stim_elecs(1:2);

%% get stim timing

params = get_stim_timing(out);
%%

trials = cell(nchs);
trials_raw = cell(nchs);
artifacts = cell(nchs);
for s_ind = good_stim_elecs
    % extract all trial waveforms
    trials_sch_raw = ...
        get_stim_trials_sch(out,s_ind);
    
    trials_sch = nan(size(trials_sch_raw));
    for t = 1:size(trials_sch_raw,2)
        
        %f=figure;
        data = squeeze(trials_sch_raw(:,t,:))'; % extract EEG data for one trial as NxT matrix
        %subplot(1,2,1); plot(data');        
        
        trials_sch(:,t,:) = PREPROCESS_CCEPS(out,data,true)';

    end
                    
    good_rec_elecs = find(~isnan(A(:,s_ind)));
    for r_ind = good_rec_elecs'
        first_tp_all_trials = squeeze(trials_sch(1,:,r_ind)); % get first time point for each trial
        if mean(isnan(first_tp_all_trials)) < mt_thresh % look for nans indicating missing trials
            trials{r_ind,s_ind} = trials_sch(:,:,r_ind);
            trials_raw{r_ind,s_ind} = trials_sch_raw(:,:,r_ind);
            artifacts{r_ind,s_ind} = stim_trials_art(out,trials_sch(:,:,r_ind),s_ind);
        else
            A(r_ind,s_ind) = nan; % if too many bad trials then remove ccep
        end
    end

end

%% return info on strongest CCEP waveforms
% vectorize network, stim idx, response idx using non-nan entries
good_ccep = ~isnan(A);
% remove entries where there is no meaningful ccep
A = A(good_ccep);
stim_idx = stim_idx(good_ccep);
response_idx = response_idx(good_ccep);

% identify strongest N1s and sort into table
[A,I] = sort(A,'descend');
stim_idx = stim_idx(I);
response_idx = response_idx(I);

info = array2table([stim_idx response_idx A],'VariableNames',{'Stim','Record','Wave'});
info.StimLabel = out.chLabels(info.Stim);
info.RecordLabel = out.chLabels(info.Record);
