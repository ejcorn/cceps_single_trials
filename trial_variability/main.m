%% download and process CCEPs data, computing N1 and N2 on average waveforms - separately run this script:
%
addpath(genpath('.'));
do_run_eli;

%% redownload CCEPs data at trial level, process with low pass filter + demeaning, compute N1 and N2 for each trial
get_stim_trials_all_pts;

compute_trial_summary_metrics;

%% assign each electrode to an anatomical parcel in a few atlases (will use Brainnetome throughout though)
parcellate;

%% compute measures of trial-by-trial variability and save for analysis in R
save_spearman_locs_soz;

%}

%% compute correlation between pre-stim phase and N1/N2
%

prestim_hippocampal_phase_vs_waveform_clcorr;

prestim_local_phase_vs_waveform_clcorr;

%}

%% plot schematics/figures for Fig 1 and fig 2

fig1_ts;
plot_ccep_waveform_fig1;
plot_ccep_waveform_fig2;