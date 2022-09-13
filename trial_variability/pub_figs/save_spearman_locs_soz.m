clear all; close all; clc
addpath(genpath('.'));
%% define patient and save directories

locations = cceps_files;
addpath(genpath(locations.freesurfer_matlab));
datadir = fullfile(locations.results_folder,'trial_metrics'); 
savedir = fullfile(locations.results_folder,'pub_figs'); mkdir(savedir);
coords = load('../data/elecs.mat');

%% Set some parameters

p_thresh = 0.05; % significance threshold
waveforms = {'N1','N2','PreStim'};
nwaves = length(waveforms);
%wave = 'N1';

%% load SOZ data
soz_data = readtable(fullfile(locations.data_folder,'SOZ.csv'));
soz_data.Properties.RowNames = soz_data.name;
soz_data = soz_data(ismember(soz_data.name,locations.subjects),:);
soz_data = soz_data(~isempty_c(soz_data.SOZElectrode),:);

%% load interictal spiking data

spike_data = load(fullfile(locations.data_folder,'spikes.mat'));
spike_quality = readtable(fullfile(locations.data_folder,'Spikes.csv')); 

%% load WM data

%wm_mask = niftiread(fullfile(locations.data_folder,'nifti','whitemattermask.nii'));
wm_mask = load_nifti(fullfile(locations.data_folder,'nifti','whitemattermask.nii'));
wm_mask.pixdim = wm_mask.pixdim*2;

%% define output structs

spear_results = struct();
elec_info = struct();
ave_cceps = struct();

%% loop through waveforms then patients, compute spearmans, save for all pts
%pt = 'HUP213';
for pt = locations.subjects
    pt = char(pt);          
    disp(pt)
    %% load data
    % trial level data
    load(fullfile(locations.results_folder,'all_trials',[pt,'_CCEPTrialWaveforms.mat']),'trials');
    load(fullfile(datadir,[pt,'_CCEPTrialSummaryMetrics.mat']));

    % load out struct and other patient data, add coordinates
    load([locations.results_folder,'/results_',pt,'_CCEP.mat']);
    clinical = pull_clinical_info([pt,'_CCEP']);    

    out = ADD_COORDINATES_BIPOLAR(out,coords);
    out.chLabels_ana = anatomic_location(out.chLabels,clinical,1);
    
    if isempty(out.chLabels_ana)
        %
        out.chLabels_ana = coords.info(find(strcmp({coords.info.name},pt))).elecs.anatomy;
        out.chLabels_ana(isempty_c(out.chLabels_ana)) = {''};

        % keeping this in here
        out.chLabels_ana = [repmat({'ECG EEG'},[8 1]);out.chLabels_ana];

        if isempty(out.chLabels_ana) | length(out.chLabels) ~= length(out.chLabels_ana)
            out.chLabels_ana = out.chLabels;       
        end
    end
    %}
    
    % get stim timing
    params = get_stim_timing(out);
    for wave = waveforms
        wave = char(wave);            

        %% quantify monotonic change with in waveform increasing trials

        [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_SEQ_SPEARMAN(out,results.(wave).data);
        spear_results.(wave).(pt).spear_trials = spear_trials;
        spear_results.(wave).(pt).spear_trials_p_adj = spear_trials_p_adj;
        
        fprintf('# significant results %s %s: %d\n',wave,pt,sum(spear_trials_p_adj<0.05,'all'))
        spear_results.(wave).(pt).good_cceps = ~isnan(spear_trials);
        
        spear_results.(wave).(pt).fail_rate = COMPUTE_FAILURE_RATE(out,results.(wave).data_0);
        [spear_results.(wave).(pt).coefvar_trials,spear_results.(wave).(pt).std_trials] = ...
            COMPUTE_CROSS_TRIAL_COEFVAR(out,results.(wave).data);
        
    end
    %%
    %{
    figure; plot(spear_results.PreStim.(pt).fail_rate(spear_results.PreStim.(pt).good_cceps),...
        spear_results.N2.(pt).fail_rate(spear_results.PreStim.(pt).good_cceps),'.');
    
    figure; plot(spear_results.PreStim.(pt).fail_rate(spear_results.PreStim.(pt).good_cceps),...
        spear_results.N1.(pt).fail_rate(spear_results.PreStim.(pt).good_cceps),'.');
    %
    high_fail = find(spear_results.N1.(pt).fail_rate>0.2);
    for j = 1:min([5 length(high_fail)])
        fail_trial_idx = find(results.N1.data_0{high_fail(j)} == 0);
        success_trial_idx = find(results.N1.data_0{high_fail(j)} > 0);
        f = figure;
        subplot(1,2,1);
        PLOT_TRIALS_WAVE_WINDOW(params,trials{high_fail(j)}(:,fail_trial_idx),params.n1_time);
        title('Fails');
        subplot(1,2,2);        
        PLOT_TRIALS_WAVE_WINDOW(params,trials{high_fail(j)}(:,success_trial_idx),params.n1_time);
        title('Success');
    end
    PreStimVN2 = cellfun(@(x,y) corr(x,y),results.PreStim.data_0(~isnan(spear_trials)),...
        results.N2.data_0(~isnan(spear_trials)));
    %figure; histogram(PreStimVN2)
    %}
    %% get bipolar array coords to see if electrodes are in WM
    
    arr = nan(size(out.locs_bipolar));
    [arr(:,1), arr(:,2), arr(:,3)] = mni2orFROMxyz(out.locs_bipolar(:,1),out.locs_bipolar(:,2),out.locs_bipolar(:,3),2,'mni');
    elec_info.(pt).arr = arr;
    elec_info.(pt).mni = out.locs_bipolar;
    elec_info.(pt).D = squareform(pdist(out.locs_bipolar,'euclidean'));

    % parcellate each electrode based on WM mask to assign each electgrodes
    % as WM or not
    [~,elec_is_wm,wm_electrodes,~,~] = ...
    CCEPS_TO_PARCELS_MNI(elec_info.(pt).mni,ones(size(elec_info.(pt).mni,1)),wm_mask,1,1);
    elec_is_wm = elec_is_wm == 1;

    % save wm electrodes
    elec_info.(pt).wm = elec_is_wm;
    
    %% save channel labels
    
    elec_info.(pt).chLabels = out.chLabels;
    elec_info.(pt).chLabels_ana = out.chLabels_ana;
    
    %% use brainnetome atlas to identify when electrodes are specifically in GM regions
        
    elec_info.(pt).gm = ~isnan(out.parcellation.elec_parcels.Brainnetome);
    elec_info.(pt).gm_tight = ~isnan(out.parcellation_tight.elec_parcels.Brainnetome);
    
    
    %% save brainnetome parcels
    
    elec_info.(pt).Brainnetome = out.parcellation.elec_parcels.Brainnetome;
    
    %% save start times
    
    elec_info.(pt).start_time = nan(length(out.stim_chs),1);
    elec_info.(pt).start_time(find(~isempty_c({out.other.periods.start_time}))) = [out.other.periods.start_time];
    [~,I] = sort(elec_info.(pt).start_time,'MissingPlacement','last');
    elec_info.(pt).chLabel_stim_order = out.chLabels(I);
    
    %% get soz electrode mask -- if SOZ data available for patient
    
    if any(strcmp(soz_data.Properties.RowNames',pt))
        soz_elec_names = soz_data{pt,'SOZElectrode'}; % get list of electrodes involved in SOZ
        soz_elec_names = strtrim(strsplit(soz_elec_names{:},',')); % split list by commas, trim white space

        % find soz electrode names in channel labels to get soz mask
        elec_info.(pt).soz = ismember(out.chLabels,soz_elec_names);             

        %% get minimum distance to any SOZ electrode

        soz_centroid = mean(out.locs_bipolar(elec_info.(pt).soz,:));
        soz_dist_tmp = nan(size(out.chLabels));
        for j = 1:length(out.chLabels)
            elec_loc = out.locs_bipolar(j,:);
            soz_ed = sqrt(sum((out.locs_bipolar(elec_info.(pt).soz,:) - elec_loc).^2,2)); % euclidean distance from electrode j to each of seizure electrodes
            soz_dist_tmp(j) = min(soz_ed); % save the closest distance to any SOZ electrode
        end

        elec_info.(pt).soz_dist = soz_dist_tmp;
    else
        elec_info.(pt).soz = nan(size(elec_is_wm));
        elec_info.(pt).soz_dist = nan(size(elec_is_wm));
    end
    
    %% save interictal spike rates for each electrode
    
    spikes_all_elecs = nan(size(out.chLabels));
    
    % only access spike data if the PPV for manual validation is good
    if spike_quality.PPV_car_(strcmp(pt,spike_quality.name)) > 0.7
        spikes = spike_data.alfredo(strcmp(pt,{spike_data.alfredo.name})).spikes;
        chLabels_spikes = spike_data.alfredo(strcmp(pt,{spike_data.alfredo.name})).labels;                
        for ch = 1:length(out.chLabels) 
            % loop through channels with CCEPs data, match channel names to
            % index in spike data, save spike data for corresponding
            % electrode in CCEP montage
            ch_idx = find(ismember(chLabels_spikes,out.chLabels(ch)));
            if ~isempty(ch_idx)
                spikes_all_elecs(ch) = spikes(ch_idx);
            end
        end
    end
    
    elec_info.(pt).spikes = spikes_all_elecs;
    
    %% save average CCEPs data
    
    ave_cceps.(pt).N1 = out.network(1).A;
    ave_cceps.(pt).N2 = out.network(2).A;
    
    % look at within SOZ CCEPs - there are hardly any
    %out.network(1).A(elec_info.(pt).soz,elec_info.(pt).soz)
    %out.network(2).A(elec_info.(pt).soz,elec_info.(pt).soz)
end

save(fullfile(savedir,'AllPatientsSpearmanLocsSoz.mat'),'spear_results','elec_info','ave_cceps');
