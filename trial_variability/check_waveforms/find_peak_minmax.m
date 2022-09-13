function [peak,peak_time] = find_peak_minmax(eeg,params,peak_range,nopeakval)

% this function will find both local minima and local maxima
% then take the maximum absolute value peak

%{
suggested values:
stim_idx = elecs(ich).stim_idx % not really a suggestion
peak_range = [50e-3 300e-3];
idx_before_stim = 30;
fs = stim.fs % not really a suggestion
%}

if ~exist('nopeakval','var')
    nopeakval=0;
end

fs = params.fs;
stim_idx = params.stim_idx;
idx_before_stim = params.idx_before_stim;
index_range = round(peak_range*fs + stim_idx);

nchs = size(eeg,2);
peak = nan(nchs,1);
peak_idx = nan(nchs,1);
for ich = 1:nchs
    values = eeg(:,ich);
    
    if all(isnan(values)) % if the trial is artifactual, it should always be NaN and never be 0
        peak(ich) = NaN;
        peak_idx(ich) = NaN;
    else
        values_in_range = values(index_range(1):index_range(2));
        baseline = mean(values(1:idx_before_stim));
        baseline_diff = (values_in_range-baseline);

        [pks_pos,locs_pos] = findpeaks(abs(baseline_diff),'MinPeakDistance',5e-3*fs); % find local maxima
        [pks_neg,locs_neg] = findpeaks(-abs(baseline_diff),'MinPeakDistance',5e-3*fs); % find local minima
        pks = [pks_pos;-pks_neg];
        locs = [locs_pos;locs_neg]; % flip sign of negative peaks again so units are positive (abs val should always be +)
        if isempty(pks)
            peak(ich) = nopeakval;
            peak_idx(ich) = NaN;
        else
            [max_pk,I] = max(pks); % find the biggest
            peak(ich) = max_pk;
            peak_idx(ich) = (round(locs(I)) + index_range(1) - 1);
        end
    end
        
    % Plot the peak
    if 0
        plot(values)
        hold on
        plot(peak_idx(ich),values(peak_idx(ich)),'o')
        plot(xlim,[baseline baseline],'r--')
        plot([stim_idx stim_idx],ylim,'r-');
        plot([index_range(1) index_range(1)],ylim,'k--')
        plot([index_range(2) index_range(2)],ylim,'k--')
        hold off
        pause
    end
    
end


peak_time = (peak_idx-stim_idx)/fs;

end