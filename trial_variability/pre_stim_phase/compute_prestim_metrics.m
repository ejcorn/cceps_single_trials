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

pt = 'HUP213';
for pt = locations.subjects
    pt = char(pt);
    disp(pt);
    load([locations.results_folder,'/results_',pt,'_CCEP.mat']);
    load(fullfile(datadir,[pt,'_CCEPTrialWaveForms.mat']),'trials','artifacts','info');
    clinical = pull_clinical_info([pt,'_CCEP']);    
    
    out = ADD_COORDINATES_BIPOLAR(out,coords);
    out.chLabels_ana = anatomic_location(out.chLabels,clinical,1);   
    
    %% load pre stim data
    
    load(fullfile(savedir,[pt,'_CCEPPrestimData.mat']));
    
    
    %% get stim params
    
    params = get_stim_timing(out);
    
    %% only use good cceps
    % require that CCEPs not excluded by earlier preprocessing,
    % and both N1 and N2 are suprathreshold.
    good_cceps = ~isempty_c(trials) & out.network(1).A>0 & out.network(2).A>0;
      
    %% compute pre-stim metrics at recording electrodes
    
    results = struct();
    %load(fullfile(savedir,[pt,'_CCEPPrestimMetrics.mat']),'results');
    
    freqBands = set_freqband_params;
    time_to_keep = ceil(params.prestim_short*out.other.stim.fs);
    %{
    for band = 1:length(freqBands.BandNames)
        prestim_metric_rch = cell(length(trials));
        bandName = freqBands.BandNames{band};
        for sch = 1:length(trials) % loop through stim electrodes
            for rch = 1:length(trials) % loops through recording electrodes
                if good_cceps(rch,sch)    
                    % compute pre-stim functional measures AT RECORDING ELECTRODE
                    data = trials{rch,sch}';
                    data = data(:,params.pre_stim_idx);
                    ft_data = fieldtrip_format(data, out.other.stim.fs,[], []);
                    Phi = GET_INST_PHASE(ft_data,freqBands.BandRanges{band});
                    prestim_metric_rch{rch,sch} = Phi(:,(end-time_to_keep):end);
                end
            end
        end
        results.(bandName).Record.data = prestim_metric_rch;
    end
    %}
    %% compute pre-stim metrics at STIM electrode
    
    % pre-populate results cell
    %{
    for band = 1:length(freqBands.BandNames)
        bandName = freqBands.BandNames{band};
        results.(bandName).Stim.data = cell(length(trials),1);
        results.(bandName).Stim.trials = cell(length(trials),1);
    end
    %}
    % loop through stim trials, download data, compute phase for each band
    all_traces_cell = cell([length(trials) 1]);
    for sch = 1:length(trials) % loop through stim electrodes
            % compute pre-stim functional measures AT STIM ELECTRODE
        if any(good_cceps(:,sch)) % if there are any good cceps for that stim electrode
            all_traces_sch = get_stim_trials_sch(out,sch,false); % download data
            all_traces_cell{sch} = all_traces_sch;
            
            % visualize whole montage time series during stimulation of one channel
            %{
            X = reshape(all_traces,[size(all_traces,1)*size(all_traces,2) size(all_traces,3)]);            
            
            f=figure;
            
            subplot(1,3,1);
            plot(zscore(X) + [1:size(X,2)])
            subplot(1,3,2);
            plot(X(:,sch));
            subplot(1,3,3);
            
            chs = [(sch-3):(sch+3)];
            plot(zscore(X(:,chs)) + [1:7]*10);
            yticklabels(out.chLabels(chs)); yticks([1:7]*10); title(sprintf('Stim Electrode %s',out.chLabels{sch}));
            xlabel('Time'); ylabel('Channel');
            f=FIGURE_SIZE_CM(f,18,9);
            saveas(f,fullfile(savedir,'ts_plot',sprintf('%sExampleTimeSeriesSurroundingStimElectrode%s.pdf',pt,out.chLabels{sch})));
            %}
            
            %{
            data = squeeze(all_traces(params.pre_stim_idx,:,sch))'; % select pre stim window
            results.(bandName).Stim.trials{sch} = data';
            if ~all(isnan(data))
                ft_data = fieldtrip_format(data, out.other.stim.fs,[], []); % convert to fieldtrip

                for band = 1:length(freqBands.BandNames) % loop thrugh bands
                    bandName = freqBands.BandNames{band};
                    Phi = GET_INST_PHASE(ft_data,freqBands.BandRanges{band}); % compute phase for each band
                    results.(bandName).Stim.data{sch} = Phi(:,(end-time_to_keep):end); % store in results struct
                end           
            end
            %}
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
    
    %% store processed data from pre
    results.data = all_traces_processed_eeg_prestim;
    
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
    %% Prep prestim data for PCA - exclude bad trials and stim electrodes
        
    % join time series from all stim electrodes into one
    X = vertcat(all_traces_reshaped_processed_eeg_prestim{~isempty_c(all_traces_reshaped_processed_eeg_prestim)});   
    X_unproc =  vertcat(all_traces_reshaped_prestim{~isempty_c(all_traces_reshaped_prestim)});
    
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
    trial_len = size(all_traces_sch_proc_prestim,1);
    n_trials = size(X,1)/trial_len;
    trial_idx = repelem(1:n_trials,trial_len)';        
    
    % exclude bad trials
    tp_nan_thresh = 0; %remove any trials with any nans in any electrode
    % find trials (trial_idx == x) with any nans in it
    bad_trials = cell2mat(cellfun(@(x) mean(isnan(X(trial_idx==x,excl_elecs)),[1 2]) > tp_nan_thresh,num2cell(1:n_trials),'UniformOutput',false))';
    % exclude the time points belonging to those trials
    excl_tps(ismember(trial_idx,find(bad_trials))) = 0;
    
    %any(~any(isnan(X(:,excl_elecs)),2))
    
    % plot distribution of missingness in rows and columns - is it bad
    % elecs or bad trials?
    
    %
    Y = X(excl_tps,excl_elecs);
    f=figure;
    subplot(1,4,1);
    imagesc(isnan(X(:,excl_elecs))); colorbar;
    title('Pre time point exclusion');
    subplot(1,4,2);
    histogram(mean(isnan(Y),1));
    title('n missing per col')
    subplot(1,4,3);
    histogram(mean(isnan(Y),2));
    title('n missing per row')    
    subplot(1,4,4);
    imagesc(isnan(X(excl_tps,excl_elecs))); colorbar;
    title('Post time point exclusion');    
    
    f=figure;
    subplot(1,2,1);
    TS_LINE_PLOT(out,X_unproc(~excl_tps,:),10,[],[],'plot all')
    subplot(1,2,2);
    TS_LINE_PLOT(out,X(excl_tps,excl_elecs),10,[],[],'plot all')
    y = X_unproc(~excl_tps,:);
    %}
    %% perform PCA on pre stim data
        
    X_fit = X(excl_tps,excl_elecs);
    [pc.coeff,pc.scores,pc.latent,pc.tsquared,pc.explained,pc.mu] = ...
    pca(X_fit);
    
    n_pcs = 20;
    f=figure; 
    subplot(1,4,1);
    imagesc(X_fit); 
    sd = std(X_fit,[],'all');
    u = mean(X_fit,'all');
    caxis([u-2*sd u+2*sd]);
    title('Data matrix');
    xlabel('channel');
    ylabel('time');

    subplot(1,4,2);
    imagesc(pc.coeff(:,1:n_pcs)');
    title('Coefficients')
    xlabel('channel');
    ylabel('component');

    subplot(1,4,3);
    legend_labels = cellfun(@(x) ['PC ',num2str(x)],num2cell(1:4)','UniformOutput',false);
    plot(pc.scores(:,1:4)); legend(legend_labels)
    xlabel('time');
    ylabel('score (a.u.)')

    subplot(1,4,4);
    plot(cumsum(pc.explained(1:n_pcs)),'b');
    ylabel('Cumulative Variance Explained');
    xlabel('Number of components');

    
    
    % you need to first exclude channels with all missing data (bipolar
    % edges)
    % then exclude time points with missing data?
    % stim electrodes probably all nans too so just get rid of those
    % electrodes
    
    % CALC TIME FREQ ANALYSIS on pre stim data - returning all nans but
    % check deltapower_windowstep_base.m in stim_predict repo
    
    %{
    d_trial = squeeze(all_traces(params.pre_stim_idx,1,:));
    [data,ft_data] = PREPROCESS_EEG(d_trial',out.other.stim.fs,out.chLabels,true);    
    ft_data.label = out.chLabels;
    res = TIMEFREQ_MULTITAPER(ft_data,freqBands.foi_hifreq,out.chLabels);
    pow = squeeze(res.powspctrm);
    
    f=figure;
    imagesc(squeeze(pow(:,:,45)))
    %}
    save(fullfile(savedir,[pt,'_CCEPPrestimMetrics.mat']),'results');
    
end
    
