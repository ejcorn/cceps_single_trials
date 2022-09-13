function all_traces = get_stim_trials_sr_pair(out,sch,rch)

%% Parameters
do_bipolar = 1;
stim_time = [-5e-3 15e-3];

%% Get various path locations
locations = cceps_files; % Need to make a file pointing to you own path
pwfile = locations.pwfile;
loginname = locations.loginname;
script_folder = locations.script_folder;

% add paths
addpath(genpath(script_folder));
if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end

%% Get info about time to download
chLabels = out.chLabels(:,1);
if ischar(sch)
    sch = find(strcmp(chLabels,sch));
end
if ischar(rch)
    rch = find(strcmp(chLabels,rch));
end
name = out.name;
fs = out.other.stim.fs;
if ~isfield(out,'clinical')
    start_time = 1;
elseif isnan(out.clinical.start_time)
    start_time = 1;
else
    start_time = out.clinical.start_time;
end

stim_indices = round(stim_time(1)*fs):round(stim_time(2)*fs);

arts = sort(out.elecs(sch).arts);

first_stim_time = arts(1)/fs;
last_stim_time= arts(end)/fs;
surround_times = out.elecs(sch).times;
download_times = [first_stim_time+surround_times(1)+start_time,last_stim_time+surround_times(2)+start_time];
arts_rel_file = arts/fs + start_time;
rel_arts = round((arts_rel_file-download_times(1))*fs);

%% Download ieeg
data = download_eeg(name,loginname, pwfile,download_times); % data object
values = data.values; % eeg time series
bits = [rel_arts + repmat(surround_times(1)*fs,length(rel_arts),1),...
    rel_arts + repmat(surround_times(2)*fs,length(rel_arts),1)];
bits = round(bits);
bits(1,1) = 1;
bits(end,2) = size(values,1);



avg = squeeze(out.elecs(sch).avg(:,rch));
nt = size(bits,1); % get number of trials

% Get alt avg
all_traces = nan(bits(1,end)-bits(1,1)+1,nt);
cmap = colormap(parula(nt));

for t = 1:nt % loop through trials
    if do_bipolar
        vals = bipolar_montage(values(bits(t,1):bits(t,end),:),rch,chLabels);
    else
        vals = values(bits(t,1):bits(t,end),rch);
        disp('a')
    end
    if length(vals) < size(all_traces,1)
        vals = [vals;repmat(vals(end),size(all_traces,1)-length(vals),1)];
        disp('b')
    end

    if length(vals) > size(all_traces,1)
        vals(end-(length(vals)-size(all_traces,1))+1:end) = [];
        disp('c')
    end
    vals = vals-mean(vals);
    
    
    %
    all_idx = 1:length(vals);
    stim_idx = out.elecs(sch).stim_idx;
    stim_indices = stim_indices + stim_idx - 1;
    non_stim_idx = all_idx;
    non_stim_idx(ismember(non_stim_idx,stim_indices)) = [];
    
    %{
    if max(abs(vals(non_stim_idx))) > 1e3
        vals = nan(size(vals));
    end
    %}
    %
    if max(abs(vals)) > 1e3
        vals = nan(size(vals));
    end
    %}
    
    all_traces(:,t) = vals; 
    
end

end