clear all; close all; clc
addpath(genpath('.'));
%% define patient and save directories

locations = cceps_files;
datadir = fullfile(locations.results_folder,'all_trials'); 
savedir = fullfile(locations.results_folder,'trial_metrics'); mkdir(savedir);
coords = load('../data/elecs.mat');

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
      
    %% compute n1 and n2 for each trial
    
    results = struct();
    waveforms = {'N1','N2','PreStim'};
    results.('N1').function = ...
        @(x) find_peak_maxabs(x,params,params.n1_time);
    results.('N2').function = ...
        @(x) find_peak_maxabs(x,params,params.n2_time);    
    results.('PreStim').function = ...
        @(x) find_peak_maxabs(x,params,params.prestim_ctrl);
    
    % use an alternate more similar to trial averaged CCEP calculation -
    % using MATLAB findpeaks function to identify local maxima rather than
    % taking max(abs(eeg))
    
    results.('N1').function_fit = ...
        @(x) find_peak_fit(x,params,params.n1_time);
    results.('N2').function_fit = ...
        @(x) find_peak_fit(x,params,params.n2_time);
    results.('PreStim').function_fit = ...
        @(x) find_peak_fit(x,params,params.prestim_ctrl);    
    
    % same as above but make it NaN if no peak identified
    
    results.('N1').function_fit_nan = ...
        @(x) find_peak_fit(x,params,params.n1_time,nan);
    results.('N2').function_fit_nan = ...
        @(x) find_peak_fit(x,params,params.n2_time,nan);
    results.('PreStim').function_fit_nan = ...
        @(x) find_peak_fit(x,params,params.prestim_ctrl,nan); 
    
    % use peak finding algorithm that counts both local minima and maxima
    % *** this is the final one
    
    results.('N1').function_minmax_nan = ...
        @(x) find_peak_minmax(x,params,params.n1_time,nan);
    results.('N2').function_minmax_nan = ...
        @(x) find_peak_minmax(x,params,params.n2_time,nan);
    results.('PreStim').function_minmax_nan = ...
        @(x) find_peak_minmax(x,params,params.prestim_ctrl,nan); 
    
    results.('N1').function_minmax = ...
        @(x) find_peak_minmax(x,params,params.n1_time);
    results.('N2').function_minmax = ...
        @(x) find_peak_minmax(x,params,params.n2_time);
    results.('PreStim').function_minmax = ...
        @(x) find_peak_minmax(x,params,params.prestim_ctrl); 
    
    % calculate area under the curve of the period instead of trying to
    % find any peaks
    
    results.('N1').function_auc = ...
        @(x) wave_quant_auc(x,params,params.n1_time);
    results.('N2').function_auc = ...
        @(x) wave_quant_auc(x,params,params.n2_time);
    results.('PreStim').function_auc = ...
        @(x) wave_quant_auc(x,params,params.prestim_ctrl);   
    
    
    results.('N1').time = params.n1_time;
    results.('N2').time = params.n2_time;
    results.('PreStim').time = params.prestim_ctrl;
    
    out.network(3).A = out.network(1).A; % add network for prestim control
    out.network(3).which = 'PreStim';
    
    for wave = waveforms % for N1 and N2
        wave = char(wave);
        trials_waves = cell(length(trials));
        trials_waves_fit = cell(length(trials));
        trials_waves_fit_nan = cell(length(trials));
        trials_waves_fit_minmax = cell(length(trials));
        trials_waves_fit_minmax_nan = cell(length(trials));
        trials_waves_auc = cell(length(trials));
        for sch = 1:length(trials) % loop through stim electrodes
            for rch = 1:length(trials) % loops through recording electrodes
                if good_cceps(rch,sch)                    
                    trials_waves{rch,sch} = results.(wave).function(trials{rch,sch});
                    trials_waves_fit{rch,sch} = results.(wave).function_fit(trials{rch,sch});
                    trials_waves_fit_nan{rch,sch} = results.(wave).function_fit_nan(trials{rch,sch});
                    trials_waves_fit_minmax{rch,sch} = results.(wave).function_minmax(trials{rch,sch});
                    trials_waves_fit_minmax_nan{rch,sch} = results.(wave).function_minmax_nan(trials{rch,sch});
                    trials_waves_auc{rch,sch} = results.(wave).function_auc(trials{rch,sch});
                end
            end
        end        
        
        results.(wave).data_maxabs = trials_waves;
        results.(wave).data_minmax_0 = trials_waves_fit_minmax;
        results.(wave).data_minmax = trials_waves_fit_minmax_nan;
        results.(wave).data = trials_waves_fit_nan;
        results.(wave).data_0 = trials_waves_fit;
        results.(wave).data_auc = trials_waves_auc;
        
        % get distance and trial-averaged waveform strength
        results.(wave).D = squareform(pdist(out.locs_bipolar,'euclidean'));
        results.(wave).A = out.network(find(strcmp(waveforms,wave))).A;

    end
    
    save(fullfile(savedir,[pt,'_CCEPTrialSummaryMetrics.mat']),'results');
    
end
    
