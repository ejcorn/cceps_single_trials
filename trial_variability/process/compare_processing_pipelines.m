clear all; close all; clc
addpath(genpath('.'));
%% define patient and save directories

locations = cceps_files;
addpath(genpath(locations.freesurfer_matlab));
savedir = fullfile(locations.results_folder,'pub_figs','figS4_filters'); mkdir(savedir);
coords = load('../data/elecs.mat');

pt = 'HUP220';

%% process data from one stim electrode

load([locations.results_folder,'/results_',pt,'_CCEP.mat']);

good_cceps = out.network(1).A >0 & out.network(2).A>0;
great_cceps = out.network(1).A >10 & out.network(2).A>0;

stim_chs_good_cceps = find(any(good_cceps,1));
trials_raw = cell(length(out.stim_chs),1);
trials_bipolar = cell(length(out.stim_chs),1);
trials_notch = cell(length(out.stim_chs),1);
trials_lowpass = cell(length(out.stim_chs),1);
trials_highpass = cell(length(out.stim_chs),1);

%% download eeg data
j = 13;
s_ind = stim_chs_good_cceps(j);
    % extract all trial waveforms
trials_sch_raw = ...
    get_stim_trials_sch(out,s_ind);
trials_raw{s_ind} = trials_sch_raw;

%% plot monopolar vs bipolar

X = squeeze(trials_sch_raw(:,20,:));
Xf = FILTER_CCEPS(X,out.other.stim.fs);

j = 100;

figure; plot(zscore(Xf(:,j+1) - Xf(:,j))); hold on; plot(zscore(Xf(:,j))); title('j+1 - j');
figure; plot(zscore(Xf(:,j) - Xf(:,j+1))); hold on; plot(zscore(Xf(:,j))); title('j - j+1');

%corr(X(:,3)-X(:,2),X(:,2)-X(:,3))
%% process EEG data

trials_bipolar_tmp = nan(size(trials_sch_raw));
trials_notch_tmp = nan(size(trials_sch_raw));
trials_lowpass_tmp = nan(size(trials_sch_raw));
trials_highpass_tmp = nan(size(trials_sch_raw));

for t = 1:size(trials_sch_raw,2)

    %f=figure;
    data = squeeze(trials_sch_raw(:,t,:))'; % extract EEG data for one trial as NxT matrix
    %subplot(1,2,1); plot(data');        
    trials_bipolar_tmp(:,t,:) = bipolar_montage(data',[],out.chLabels);
    trials_lowpass_tmp(:,t,:) = PREPROCESS_CCEPS(out,data,true)';
    trials_highpass_tmp(:,t,:) = PREPROCESS_CCEPS_HIGHPASS(out,data,true)';
    trials_notch_tmp(:,t,:) = PREPROCESS_CCEPS_LINE_ONLY(out,data,true)';

end

trials_bipolar{s_ind} = trials_bipolar_tmp;
trials_notch{s_ind} = trials_notch_tmp;
trials_lowpass{s_ind} = trials_lowpass_tmp;
trials_highpass{s_ind} = trials_highpass_tmp;

%% plot

% specify which stim channel, how many CCEPs, which trial
j = 1;
n_cceps_plot = 10; % how many CCEPs to plot
trial_plot = 15; % which trial to plot (shouldn't matter

s_ind = stim_chs_good_cceps(j);
r_inds = find(great_cceps(:,s_ind));
r_inds = r_inds(1:min(length(r_inds),n_cceps_plot));

% extract data from different preprocessing pipelines
data_bipolar = squeeze(trials_bipolar{s_ind}(:,trial_plot,r_inds));
data_notch = squeeze(trials_notch{s_ind}(:,trial_plot,r_inds));
data_lowpass = squeeze(trials_lowpass{s_ind}(:,trial_plot,r_inds));

% stack and scale the data so it can be plotted on one axis
data_mean = mean(data_lowpass);
data_sd = std(data_lowpass);
sf = 10;

offset_matrix = sf*repmat(1:size(data_lowpass,2),[size(data_lowpass,1) 1]);
data_bipolar = data_bipolar./data_sd + offset_matrix;
data_notch = data_notch./data_sd  + offset_matrix;
data_lowpass = data_lowpass./data_sd  + offset_matrix;

% make time matrix
t_plot = [-200 400]; % -200 to 400 ms relative to stim
times = linspace(-500,800,size(data_bipolar,1))';
times = repmat(times,[1 length(r_inds)]);

% make space for n1 and n2 windows
window1 = params.n1_time;
window2 = params.n2_time;
rec1 = @(x) rectangle('Position',[1000*window1(1) -10000 1000*diff(window1) 50000],...
    'FaceColor',[0.3 0.3 0.3 0.3],'EdgeColor',[0.3 0.3 0.3 0]);
rec2 = @(x) rectangle('Position',[1000*window2(1) -10000 1000*diff(window2) 50000],...
    'FaceColor',[0.9 0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9 0]);

% compute y limits

data_all = [data_bipolar data_notch data_lowpass];
yl = [min(data_all,[],[1 2]) max(data_all,[],[1 2])];
yl = [yl(1) - 0.1*abs(diff(yl)), yl(2) + 0.1*abs(diff(yl))];

% make plot

f = figure;
subplot(1,3,1);
rec1(); hold on;
rec2();
plot(times,data_bipolar);
xlim(t_plot)
ylim(yl)
title('Unfiltered');
yticks(offset_matrix(1,:)); yticklabels(out.chLabels(r_inds));
xlabel('Time (ms)');
prettifyEJC;

subplot(1,3,2);
rec1(); hold on;
rec2();
plot(times,data_notch);
xlim(t_plot)
ylim(yl)
title('Notch Filter + Detrend')
yticks(offset_matrix(1,:)); yticklabels(out.chLabels(r_inds));
xlabel('Time (ms)');
prettifyEJC;

subplot(1,3,3);
rec1(); hold on;
rec2();
plot(times,data_lowpass);
xlim(t_plot)
ylim(yl)
title('30 Hz Low Pass + Detrend')
yticks(offset_matrix(1,:)); yticklabels(out.chLabels(r_inds));
xlabel('Time (ms)');
prettifyEJC;

f=FIGURE_SIZE_CM(f,21,9);
saveas(f,fullfile(savedir,[pt,'Stim',out.chLabels{s_ind},'Trial',num2str(trial_plot),'_CompareFilters_NoneNotchLowPass30.pdf']));

%% plot high pass

% specify which stim channel, how many CCEPs, which trial
j = 13;
n_cceps_plot = 10; % how many CCEPs to plot
trial_plot = 15; % which trial to plot (shouldn't matter

s_ind = stim_chs_good_cceps(j);
r_inds = find(great_cceps(:,s_ind));
r_inds = r_inds(1:min(length(r_inds),n_cceps_plot));

% extract data from different preprocessing pipelines
data_highpass = squeeze(trials_highpass{s_ind}(:,trial_plot,r_inds));
data_notch = squeeze(trials_notch{s_ind}(:,trial_plot,r_inds));
data_lowpass = squeeze(trials_lowpass{s_ind}(:,trial_plot,r_inds));

% stack and scale the data so it can be plotted on one axis
data_mean = mean(data_lowpass);
data_sd = std(data_lowpass);
sf = 10;

offset_matrix = sf*repmat(1:size(data_lowpass,2),[size(data_lowpass,1) 1]);
data_highpass = data_highpass./data_sd + offset_matrix;
data_notch = data_notch./data_sd  + offset_matrix;
data_lowpass = data_lowpass./data_sd  + offset_matrix;

% make time matrix
t_plot = [-200 400]; % -200 to 400 ms relative to stim
times = linspace(-500,800,size(data_highpass,1))';
times = repmat(times,[1 length(r_inds)]);

% make space for n1 and n2 windows
window1 = params.n1_time;
window2 = params.n2_time;
rec1 = @(x) rectangle('Position',[1000*window1(1) -10000 1000*diff(window1) 50000],...
    'FaceColor',[0.3 0.3 0.3 0.3],'EdgeColor',[0.3 0.3 0.3 0]);
rec2 = @(x) rectangle('Position',[1000*window2(1) -10000 1000*diff(window2) 50000],...
    'FaceColor',[0.9 0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9 0]);

% compute y limits

data_all = [data_highpass data_notch data_lowpass];
yl = [min(data_all,[],[1 2]) max(data_all,[],[1 2])];
yl = [yl(1) - 0.1*abs(diff(yl)), yl(2) + 0.1*abs(diff(yl))];

% make plot

f = figure;
subplot(1,3,1);
rec1(); hold on;
rec2();
plot(times,data_highpass);
xlim(t_plot)
ylim(yl)
title('Notch Filter + 1 Hz High Pass');
yticks(offset_matrix(1,:)); yticklabels(out.chLabels(r_inds));
xlabel('Time (ms)');
prettifyEJC;

subplot(1,3,2);
rec1(); hold on;
rec2();
plot(times,data_notch);
xlim(t_plot)
ylim(yl)
title('Notch Filter + Detrend')
yticks(offset_matrix(1,:)); yticklabels(out.chLabels(r_inds));
xlabel('Time (ms)');
prettifyEJC;

subplot(1,3,3);
rec1(); hold on;
rec2();
plot(times,data_lowpass);
xlim(t_plot)
ylim(yl)
title('30 Hz Low Pass + Detrend')
yticks(offset_matrix(1,:)); yticklabels(out.chLabels(r_inds));
xlabel('Time (ms)');
prettifyEJC;

f=FIGURE_SIZE_CM(f,21,9);
saveas(f,fullfile(savedir,[pt,'Stim',out.chLabels{s_ind},'Trial',num2str(trial_plot),'_CompareFilters_HighPassNotchLowPass30.pdf']));
close(f);