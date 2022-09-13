function [data,ft_data,data_mean] = FILTER_EEG(data,srate,elec_labels,ft,data_mean)
% INPUTS:
% data: NxT matrix of ECoG signal
% srate: sampling rate in Hertz
% elec_labels: cell array of electrode labels
% rmv: binary vector highlighting noisy electrodes with 1's... can do
% custom or calculate automatically based on data
% ft: logical offering option to convert to fieldtrip format
% data_demean: vector of input. default is to use input data

%{
if ~exist('rmv','var')
    rmv = reject_elecs(data, 2); % thr of 2 stds, %data-driven: find channels with large kurtosis or linelength
elseif exist('rmv','var')
    if isempty(rmv) % if you give blank input
        rmv = reject_elecs(data, 2); % thr of 2 stds, %data-driven: find channels with large kurtosis or linelength
    end
end

%data = data(~rmv,:);
%elec_labels = elec_labels(~rmv,:);
%}
%% Look at line noise

% filter out 60 Hz harmonics
% first input is order, then bounds of filter, and type of filter

ch_msg = ~any(isnan(data),2);

%
[b,a] = butter(4, [56.5/(srate/2), 63.5/(srate/2)], 'stop');
data(ch_msg,:) = filtfilt(b,a,data(ch_msg,:)')'; % filtfilt prevents phase distortion

[b,a] = butter(4, [116.5/(srate/2), 123.5/(srate/2)], 'stop');
data(ch_msg,:) = filtfilt(b,a,data(ch_msg,:)')';

[b,a] = butter(4, [176.5/(srate/2), 183.5/(srate/2)], 'stop');
data(ch_msg,:) = filtfilt(b,a,data(ch_msg,:)')';
%}

%% demean - skip for now and do after bipolar referencing
%{
if  ~exist('data_mean','var')
	data_mean = mean(data,2);
elseif isnan(data_mean)
    data_mean = repmat(0,[size(data,1) 1]); % don't demean if set to nan
end

data = data - repmat(data_mean,[1 size(data,2)]); % demean
%}
%% detrend

%data = detrend(data')'; % as opposed to high pass filtering

%% CAR - common average reference
% if data looks to have different levels of noise based on elec, might want
% to do this in groups
%
% only perform CAR on electrodes
%{
data(~rmv,:) = get_CAR(data(~rmv,:), elec_labels(~rmv,:)); % CAR by group - it expects elec labels to be strings, where the letters define the electrode and the number defines the contact
data(rmv,:) = NaN; % replace data from removed electrodes with NaNs
%}
%% Get into fieldtrip format
% last argument is for events structure
if ~exist('ft','var') % default 
    ft = false;
end
ft_data = []; % by default don't convert data, return empty variable
if ft
    ft_data = fieldtrip_format(data, srate, elec_labels, []);
end
%}