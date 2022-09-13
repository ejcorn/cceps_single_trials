%% this script downloads data and processes it for pre stim spectral analysis
clear all; close all; clc
addpath(genpath('.'));
%% define patient and save directories

locations = cceps_files;
addpath(locations.fieldtrip); 
ft_defaults;
datadir = fullfile(locations.results_folder,'all_trials'); 
savedir = fullfile(locations.results_folder,'prestim_metrics'); mkdir(savedir);
mkdir(fullfile(savedir,'ts_plot'));
coords = load('../data/elecs.mat');

pt = 'HUP212';
for pt = locations.subjects
    pt = char(pt);
    disp(pt);
    load([locations.results_folder,'/results_',pt,'_CCEP.mat']);
    load(fullfile(datadir,[pt,'_CCEPTrialWaveForms.mat']),'trials','artifacts','info');
    clinical = pull_clinical_info([pt,'_CCEP']);    

    out = ADD_COORDINATES_BIPOLAR(out,coords);
    out.chLabels_ana = anatomic_location(out.chLabels,clinical,1);   
    
    %% get stim params
    
    params = get_stim_timing(out);
    
    %% only use good cceps
    % require that CCEPs not excluded by earlier preprocessing,
    % and both N1 and N2 are suprathreshold.
    good_cceps = ~isempty_c(trials) & out.network(1).A>0 & out.network(2).A>0;
      
    %% loop through stim trials, download data
    
    data = struct();
    all_traces_cell = cell([length(trials) 1]);    
    for sch = 1:length(trials) % loop through stim electrodes
            % compute pre-stim functional measures AT STIM ELECTRODE
        if any(good_cceps(:,sch)) % if there are any good cceps for that stim electrode
            all_traces_sch = get_stim_trials_sch(out,sch,false); % download data
            all_traces_cell{sch} = all_traces_sch;            
        end
    end
    
    %% process all pre-stim data using pipeline designed for spectral analysis
        
    % preprocess prestim eeg data then
    % reshape prestim trials to vertically concatenate them
    % so it's time x channel
    
    all_traces_reshaped = cell(size(all_traces_cell)); % store reshaped but unprocessed trial data   
    all_traces_reshaped_processed_eeg = cell(size(all_traces_cell)); % store processed and reshaped trial data   
    all_traces_reshaped_processed_cceps = cell(size(all_traces_cell)); % store processed and reshaped trial data -- using CCEPs pipeline
    
    all_traces_processed_eeg_prestim = cell(size(all_traces_cell)); % store processed but unreshaped data
    all_traces_reshaped_processed_eeg_prestim = cell(size(all_traces_cell)); % store processed and reshaped trial data -- pre stim only
    all_traces_reshaped_prestim = cell(size(all_traces_cell)); % store reshaped but unprocessed trial data -- pre stim only   
    for sch = 1:length(trials) % loop through stim electrodes
            % compute pre-stim functional measures AT STIM ELECTRODE
        if any(good_cceps(:,sch)) % if there are any good cceps for that stim electrode           
            all_traces_sch = all_traces_cell{sch};            
            all_traces_sch_proc = nan(size(all_traces_sch));
            all_traces_sch_proc_cceps = nan(size(all_traces_sch));
            for trial = 1:size(all_traces_sch,2)
                one_trial = squeeze(all_traces_sch(:,trial,:));
                one_trial_proc = PREPROCESS_EEG(out,one_trial',false)'; % process each trial separately
                all_traces_sch_proc(:,trial,:) = one_trial_proc; % store processed data in new variable
                
                all_traces_sch_proc_cceps(:,trial,:) = PREPROCESS_CCEPS(out,one_trial',false)';
            end
                        
            % ultimately not using any of these 3 for PCA - just
            % visualization as quality control
            all_traces_reshaped_processed_eeg{sch} = reshape(all_traces_sch_proc,[size(all_traces_sch_proc,1)*size(all_traces_sch_proc,2) size(all_traces_sch_proc,3)]);
            all_traces_reshaped{sch} = reshape(all_traces_sch,[size(all_traces_sch,1)*size(all_traces_sch,2) size(all_traces_sch,3)]);
            all_traces_reshaped_processed_cceps{sch} = reshape(all_traces_sch_proc_cceps,[size(all_traces_sch_proc_cceps,1)*size(all_traces_sch_proc_cceps,2) size(all_traces_sch_proc_cceps,3)]);
            
            % *** either use params.pre_stim_idx OR 1:params.idx_before_stim
            all_traces_sch_proc_prestim = all_traces_sch_proc(params.pre_stim_idx,:,:);            
            all_traces_processed_eeg_prestim{sch} = all_traces_sch_proc_prestim;
            %all_traces_sch_proc_prestim = all_traces_sch_proc(1:params.idx_before_stim,:,:);            
            all_traces_reshaped_processed_eeg_prestim{sch} = reshape(all_traces_sch_proc_prestim,[size(all_traces_sch_proc_prestim,1)*size(all_traces_sch_proc_prestim,2) size(all_traces_sch_proc_prestim,3)]);
            
            % store pre stim data only for unprocessed data - for
            % understanding why some electrodes on some trials were excluded
            all_traces_sch_prestim = all_traces_sch(params.pre_stim_idx,:,:);
            all_traces_reshaped_prestim{sch} = reshape(all_traces_sch_prestim,[size(all_traces_sch_prestim,1)*size(all_traces_sch_prestim,2) size(all_traces_sch_prestim,3)]);
            
        end
    end
    
    %% make sure preprocessing worked -- these numbers for HUP213
    % you can see 
    %
    %{
    sch = 11;
    sch_plot = 105;
    k = 5;
    
    f=figure;    
    subplot(1,4,1);
    TS_LINE_PLOT(out,all_traces_reshaped{sch},sch,k,sch_plot)
    a = gca;    
    title(sprintf('%s - Unprocessed',a.Title.String));
    subplot(1,4,2);
    TS_LINE_PLOT(out,all_traces_reshaped_processed_eeg{sch},sch,k,sch_plot)
    a = gca;    
    title(sprintf('%s - Processed EEG',a.Title.String));
    subplot(1,4,3);
    TS_LINE_PLOT(out,all_traces_reshaped_processed_eeg_prestim{sch},sch,k,sch_plot)
    a = gca;    
    title(sprintf('%s - Processed EEG - Prestim',a.Title.String));
    
    subplot(1,4,4);
    TS_LINE_PLOT(out,all_traces_reshaped_processed_cceps{sch},sch,k,sch_plot)
    a = gca;    
    title(sprintf('%s - Processed CCEPs',a.Title.String));
    %}
    
    %% store processed data from pre stim period
    
    data.raw = all_traces_cell;
    data.processed = all_traces_processed_eeg_prestim;
    data.reshape_fxn = @(y) reshape(y,[size(y,1)*size(y,2) size(y,3)]);
    
    save(fullfile(savedir,[pt,'_CCEPPrestimData.mat']),'data','-v7.3');
    
end