function data_proc = PREPROCESS_EEG(out,data,demean)

% INPUTS:
% out: cceps out struct
% data: NxT matrix of eeg data from single CCEP stimulation trial
% demean: true false - demean data relative to baseline mean?
%
% OUTPUTS:
% data: input data with following steps applied:
% - filtered pre-stim and post-stim data separately to remove line noise
% - bipolar montage
% - subtract mean of -500 ms to -60 ms relative to stim as baseline

params = get_stim_timing(out);

%rmv = reject_elecs(data(:,params.pre_stim_idx),2);
%data(rmv,:) = nan;
data = TRIALS_ARTIFACT_REJECT(data,params);
data_proc = nan(size(data));

%{
f=figure;
subplot(1,2,1);
plot(data');
%}
if ~all(isnan(data))
    % get pre and post-stim data

    prestim_data = data(:,1:(params.stim_indices(1)-1));
    intrastim_data = data(:,params.stim_indices);
    poststim_data = data(:,(1+params.stim_indices(end)):end);

    % preprocess pre-stim and post-stim data - filter line noise 

    prestim_data = FILTER_EEG(prestim_data,out.other.stim.fs); 
    poststim_data = FILTER_EEG(poststim_data,out.other.stim.fs);%,[],[],data_mean);

    % store processed data in matrix

    if size([prestim_data intrastim_data poststim_data],2) ~= size(data,2)
        error('DIMENSION MISMATCH')
    end

    data_proc(:,1:(params.stim_indices(1)-1)) = prestim_data;
    data_proc(:,params.stim_indices) = intrastim_data;
    data_proc(:,(1+params.stim_indices(end)):end) = poststim_data;  

    % do bipolar montage                    
    data_proc = bipolar_montage(data_proc',[],out.chLabels)';

    % demean bipolar montage data relative to prestim bipolar data
    if demean
        prestim_bipolar = bipolar_montage(prestim_data(:,1:params.idx_before_stim)',[],out.chLabels);
        prestim_mean = mean(prestim_bipolar,1);
        data_proc = data_proc - repmat(prestim_mean',[1 size(data_proc,2)]);
    end
    
    %subplot(1,2,2); plot(data_proc');
    
    %{
    % look at correlation structure pre and post stim for each trial
    f=figure; 
    subplot(1,2,1); imagesc(corr(prestim_data')); caxis([-1 1]);
    subplot(1,2,2); imagesc(corr(poststim_data')); caxis([-1 1]);
    %}   
else
    % if all data are nans just put it in the all trial matrix
    data_proc = data;
end