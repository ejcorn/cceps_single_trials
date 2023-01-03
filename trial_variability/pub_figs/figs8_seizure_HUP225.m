clear all; close all; clc
addpath(genpath('.'));
%% define patient and save directories

locations = cceps_files;
addpath(genpath(locations.freesurfer_matlab));
savedir = fullfile(locations.results_folder,'pub_figs','figS8_seizures'); mkdir(savedir);
coords = load('../data/elecs.mat');

pt = 'HUP225';

%% load SOZ data
soz_data = readtable(fullfile(locations.data_folder,'SOZ.csv'));
soz_data.Properties.RowNames = soz_data.name;
soz_data = soz_data(ismember(soz_data.name,locations.subjects),:);
soz_data = soz_data(~isempty_c(soz_data.SOZElectrode),:);
soz_pt = soz_data{pt,'SOZElectrode'};

%% process data from one stim electrode

load([locations.results_folder,'/results_',pt,'_CCEP.mat']);
load(fullfile(locations.results_folder,'pub_figs','AllPatientsSpearmanLocsSoz.mat'));
load(fullfile(locations.results_folder,'pub_figs','fig4','AllPatientsHippocampalPhaseCLcorr.mat'))
load(fullfile(locations.results_folder,'all_trials',[pt,'_CCEPTrialWaveforms.mat']),'trials');

seizure_elec = 'RC1';
seizure_elec_ind = find(strcmp(seizure_elec,out.chLabels));

%% are monotonic trends in HUP 225 occurring during seizure?

sig_monotonic_all = struct();
for wave = {'N1','N2'}
    wave = char(wave);
    disp(wave);
    disp(sum(spear_results.(wave).(pt).spear_trials_p_adj(:,seizure_elec_ind)<0.05))
    sig_monotonic = spear_results.(wave).(pt).spear_trials_p_adj(:,seizure_elec_ind)<0.05;
    disp(spear_results.(wave).(pt).spear_trials(sig_monotonic,seizure_elec_ind))
    disp(out.chLabels(sig_monotonic))
    sig_monotonic_all.(wave) = sig_monotonic;
end

%% plot trials with significant n1 trends preceding seizure

% extract trials into time series of trials in order for each channel
params = get_stim_timing(out);
r_ind = find(sig_monotonic_all.N1);
ts = cell(length(r_ind),1);

trials_plot = 1:1:30;
for j = 1:length(r_ind)
    X = trials{r_ind(j),seizure_elec_ind};
    % extract -500 to 500 ms epochs for each stim
    % throw out every other trial to make visualization easier
    X = X(params.trial_even_plot_idx,trials_plot);
    ts{j} = reshape(X,[],1);
end
    
ts = [ts{:}];

% plot trials

f=figure;
TS_LINE_PLOT_SIMPLE(out,ts,seizure_elec_ind,r_ind)
title([seizure_elec,' Stim Preceding Focal Aware Seizure']);
prettifyEJC;
f=FIGURE_SIZE_CM(f,21,9);
saveas(f,fullfile(savedir,[pt,'PreSeizureStimTrials30.pdf']));
close(f);

% repeat above plot but only show 1/3 of trials
ts = cell(length(r_ind),1);

trials_plot = 1:3:30;
for j = 1:length(r_ind)
    X = trials{r_ind(j),seizure_elec_ind};
    % extract -500 to 500 ms epochs for each stim
    % throw out every other trial to make visualization easier
    X = X(params.trial_even_plot_idx,trials_plot);
    ts{j} = reshape(X,[],1);
end
    
ts = [ts{:}];

% plot trials
f=figure;
TS_LINE_PLOT_SIMPLE(out,ts,seizure_elec_ind,r_ind)
title([seizure_elec,' Stim Preceding Focal Aware Seizure']);
prettifyEJC;
f=FIGURE_SIZE_CM(f,21,9);
f.Renderer = 'painters';
print(fullfile(savedir,[pt,'PreSeizureStimTrials10.pdf']),'-dpdf');
saveas(f,fullfile(savedir,[pt,'PreSeizureStimTrials10.pdf']));
close(f);


%% is hippocampal phase locking in HUP 225 occurring during seizure?
for wave = {'N1','N2'}
    wave = char(wave);
    elecs = fieldnames(clcorr_results.(wave).(pt));

    for elec_j = 1:length(elecs)
        elec_name = elecs{elec_j};
        clcorr_r = clcorr_results.(wave).(pt).(elec_name).clcorr_trials;
        clcorr_p_adj = clcorr_results.(wave).(pt).(elec_name).clcorr_trials_p_adj;
        disp(sum(clcorr_p_adj(:,seizure_elec_ind)<0.05))

    end
end

%% 