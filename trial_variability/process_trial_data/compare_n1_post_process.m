clear all; close all; clc
addpath(genpath('.'));
%% define patient and save directories

locations = cceps_files;
addpath(genpath(locations.freesurfer_matlab));
datadir = fullfile(locations.results_folder,'trial_metrics'); 
coords = load('../data/elecs.mat');

%% loop through waveforms then patients, compute spearmans, save for all pts
pt = 'HUP213';
results = struct();
for pt = locations.subjects
    pt = char(pt);          
    disp(pt)
    %% load data
    % trial level data
    load(fullfile(locations.results_folder,'all_trials',[pt,'_CCEPTrialWaveforms.mat']),'trials');
    % load out struct and other patient data, add coordinates
    load([locations.results_folder,'/results_',pt,'_CCEP.mat']);
    
    params = get_stim_timing(out);
    
    %% only use good cceps
    % require that CCEPs not excluded by earlier preprocessing,
    % and both N1 and N2 are suprathreshold.
    good_cceps = ~isempty_c(trials) & out.network(1).A>0 & out.network(2).A>0;
          
    %% loop through trials and compute N1 on averaged data
    results.('N1').function_fit_nan = ...
        @(x) find_peak_fit(x,params,params.n1_time,nan);
    results.('N2').function_fit_nan = ...
        @(x) find_peak_fit(x,params,params.n2_time,nan);

    results.('N1').NX_nofilter = out.network(1).A;
    results.('N2').NX_nofilter = out.network(2).A;
    for wave = {'N1','N2'}
        wave = char(wave);
        trial_avg_fit_nan = nan(length(trials));
        for sch = 1:length(trials) % loop through stim electrodes
            for rch = 1:length(trials) % loops through recording electrodes
                if good_cceps(rch,sch)                    
                    % average individual *processed* trials (i.e. they have
                    % all been filtered low pass 30 hz)
                    % to generate an average waveform
                    
                    trial_avg_processed = nanmean(trials{rch,sch},2); 
                    trial_avg_fit_nan(rch,sch) = results.(wave).function_fit_nan(trial_avg_processed);
                end
            end
        end    
        results.(wave).data = trial_avg_fit_nan;
        
    end
     
            
    f=figure('Name',pt);
    subplot(1,2,1);
    scatter(results.('N1').NX_nofilter(:),results.('N1').data(:),3,'b','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    [r p] = corr(results.('N1').NX_nofilter(:),results.('N1').data(:),'rows','pairwise');
    prettifyEJC;
    title(sprintf('%s: r = %s','N1',LABELROUND2(r)));

    subplot(1,2,2);
    scatter(results.('N2').NX_nofilter(:),results.('N2').data(:),3,'b','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    [r p] = corr(results.('N2').NX_nofilter(:),results.('N2').data(:),'rows','pairwise');
    prettifyEJC;
    title(sprintf('%s: r = %s','N2',LABELROUND2(r)));
    
    % results with 30 hz filter; most r's between 0.1 and 0.4, except
    % HUP213 was r = 0.7


end


