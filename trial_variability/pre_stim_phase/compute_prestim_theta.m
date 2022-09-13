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
    
    out = ADD_COORDINATES_BIPOLAR(out,coords);
    out.chLabels_ana = coords.info(find(strcmp({coords.info.name},pt))).elecs.anatomy;   
    
    %% load pre stim data
    
    load(fullfile(savedir,[pt,'_CCEPPrestimData.mat']));
    
    
    %% get stim params
    
    params = get_stim_timing(out);
    
    %% only use good cceps
    % require that CCEPs not excluded by earlier preprocessing,
    % and both N1 and N2 are suprathreshold.
    good_cceps = ~isempty_c(trials) & out.network(1).A>0 & out.network(2).A>0;
      
    %% compute pre-stim metrics at hippocampus
    
    results = struct();
    %load(fullfile(savedir,[pt,'_CCEPPrestimMetrics.mat']),'results');
    
    freqBands = set_freqband_params;
    time_to_keep = ceil(params.prestim_short*out.other.stim.fs);
    
    % get indices of hipocampal channels
    hippocampalChannelIdx = find(contains(out.chLabels_ana,'Hippocampus'));
    
    %
    for band = 2:4
        prestim_metric_rch = cell(length(trials),1);
        bandName = freqBands.BandNames{band};
        for sch = 1:length(trials) % loop through stim electrodes                        
            if any(good_cceps(:,sch))
                
                % compute pre-stim functional measures AT HIPPOCAMPAL RECORDING ELECTRODES
                X = data.raw{sch}(params.pre_stim_idx,:,:);                
                prestim_metric_rch{sch} = nan([size(X,1),size(X,3)]);
                
                %{
                X = squeeze(X(:,1,:));
                for j = hippocampalChannelIdx'
                    FFT_PLOT(X(:,j),out.other.stim.fs)
                end
                %}
                for t = 1:size(X,2)
                    for rch = hippocampalChannelIdx'                
                        ft_data = fieldtrip_format(squeeze(X(:,t,rch))', out.other.stim.fs,[], []);
                        Phi = GET_INST_PHASE(ft_data,freqBands.BandRanges{band});
                        prestim_metric_rch{sch}(:,rch) = Phi;
                    end
                end
                
            end
        end
        results.(bandName).Record.data = prestim_metric_rch;
    end    
    
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
    
