function bad_trials = stim_trials_art(out,all_trials,sch)

% INPUT:
% all_trials: time x number of trials matrix of ccep waveforms obtained from single trial
% out: output structure
% sch: stim channel

nt = size(all_trials,2);
bad_trials = false(nt,1); % return binary variable highlighting bad trials
for t = 1:nt % loop through trials
    
    vals = all_trials(:,nt);
    % mark trials as artifacts
    if mean(isnan(vals)) < 1 % if there is data there
        bad_trials(t) = identify_waveform_artifact(vals,sch,out);
    end
    
end

end