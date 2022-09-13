clear all; close all; clc
addpath(genpath('.'));
%% define patient and save directories

locations = cceps_files;
savedir = fullfile(locations.results_folder,'all_trials'); mkdir(savedir);

for pt = locations.subjects
    pt = char(pt);
    fprintf('Downloading trial data for %s\n',pt)
    load([locations.results_folder,'/results_',pt,'_CCEP.mat']);
    %% get single trial cceps
    opts = struct('mt_thresh',0.33);
    [trials,artifacts,info,trials_raw] = get_stim_trials_pt(out,opts);    

    %% save
    save(fullfile(savedir,[pt,'_CCEPTrialWaveForms.mat']),'trials','artifacts','info','-v7.3','-nocompression');
    % ^- note that artifacts variable is not utilized, but it contains an indicator of whether each trial
    % would meet the same artifact rejection criteria applied to the average waveforms

end