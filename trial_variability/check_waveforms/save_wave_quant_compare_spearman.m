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

%% define output structs

spear_results_wave_method = struct();

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

    out = ADD_COORDINATES_BIPOLAR(out,coords);
    
    % get stim timing
    params = get_stim_timing(out);
    for wave = waveforms
        wave = char(wave);            

        %% quantify monotonic change with in waveform increasing trials

        [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_SEQ_SPEARMAN(out,results.(wave).data_fit_0);
        spear_results_wave_method.(wave).findpeaks.(pt).spear_trials = spear_trials;
        spear_results_wave_method.(wave).findpeaks.(pt).spear_trials_p_adj = spear_trials_p_adj;
        
        fprintf('Find Peaks: # significant results %s %s: %d\n',wave,pt,sum(spear_trials_p_adj<0.05,'all'))
        
        [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_SEQ_SPEARMAN(out,results.(wave).data_fit);
        spear_results_wave_method.(wave).findpeaks_nan.(pt).spear_trials = spear_trials;
        spear_results_wave_method.(wave).findpeaks_nan.(pt).spear_trials_p_adj = spear_trials_p_adj;
        
        [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_SEQ_SPEARMAN(out,results.(wave).data_0);
        spear_results_wave_method.(wave).findpeaks_minmax.(pt).spear_trials = spear_trials;
        spear_results_wave_method.(wave).findpeaks_minmax.(pt).spear_trials_p_adj = spear_trials_p_adj;
        
        fprintf('Find Peaks MinMax: # significant results %s %s: %d\n',wave,pt,sum(spear_trials_p_adj<0.05,'all'))
        
        [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_SEQ_SPEARMAN(out,results.(wave).data);
        spear_results_wave_method.(wave).findpeaks_minmax_nan.(pt).spear_trials = spear_trials;
        spear_results_wave_method.(wave).findpeaks_minmax_nan.(pt).spear_trials_p_adj = spear_trials_p_adj;
        
        fprintf('Find Peaks MinMax NaN: # significant results %s %s: %d\n',wave,pt,sum(spear_trials_p_adj<0.05,'all'))
        
        [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_SEQ_SPEARMAN(out,results.(wave).data_auc);
        spear_results_wave_method.(wave).auc.(pt).spear_trials = spear_trials;
        spear_results_wave_method.(wave).auc.(pt).spear_trials_p_adj = spear_trials_p_adj;
        
        fprintf('AUC: # significant results %s %s: %d\n',wave,pt,sum(spear_trials_p_adj<0.05,'all'))
        
        [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_SEQ_SPEARMAN(out,results.(wave).data_maxabs);
        spear_results_wave_method.(wave).maxabs.(pt).spear_trials = spear_trials;
        spear_results_wave_method.(wave).maxabs.(pt).spear_trials_p_adj = spear_trials_p_adj;
        
        fprintf('Max Abs: # significant results %s %s: %d\n',wave,pt,sum(spear_trials_p_adj<0.05,'all'))

        
    end
    
end
   
save(fullfile(savedir,'AllPatientsWaveQuantMethodsSpearman.mat'),'spear_results_wave_method');
