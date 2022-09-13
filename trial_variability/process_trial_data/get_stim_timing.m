function params = get_stim_timing(out)

% Parameters - define indices of timing relative to stim    
    params.time_to_take = [-500e-3 800e-3]; % this is the time of each trial/epoch relative to stim
    params.n1_time = [15e-3 50e-3];
    params.loose_n1_time = [15e-3 50e-3];
    params.n2_time = [50e-3 300e-3];
    params.stim_time = [-5e-3 15e-3];
    params.tight_stim_time = [-5e-3 10e-3];
    params.stim_val_thresh = 1e3;
    params.rel_thresh = 3;
    params.fs = out.other.stim.fs;
    params.max_crossings = 3;
    
    params.prestim_ctrl = [-300e-3 -200e-3]; % control period to find peaks during pre stim period
    
    params.prestim_short = [50e-3]; % only look at prestim metrics 50 ms before stim

    params.n1_idx = floor(params.n1_time*params.fs);
    params.loose_n1_idx = floor(params.loose_n1_time*params.fs);
    params.n2_idx = floor(params.n2_time*params.fs);
    params.tight_stim_indices = floor(params.tight_stim_time*params.fs);               
    
    params.stim_indices_0 = round(params.stim_time(1)*params.fs):round(params.stim_time(2)*params.fs);
    params.all_idx = 1:size([out.elecs.avg],1); % get length of whole epoch
    params.stim_idx = unique([out.elecs.stim_idx]);    
    
    params.stim_indices = params.stim_indices_0 + params.stim_idx - 1;
    params.non_stim_idx = params.all_idx;
    params.non_stim_idx(ismember(params.non_stim_idx,params.stim_indices)) = [];
    params.pre_stim_idx = 1:min(params.stim_indices);    
    params.post_stim_idx = max(params.stim_indices):max(params.all_idx);    
    
    % to be conservative, ignore an extra 30 ms before start of stim artifact
    params.idx_before_stim = max(params.pre_stim_idx) - round(30e-3 * params.fs);

    params.trial_plot_times = [-200e-3 400e-3]; % period to plot when visualizing trial N1/N2s
    params.trial_plot_idx = params.stim_idx + floor([params.trial_plot_times(1)*params.fs : params.trial_plot_times(2)*params.fs]);

    % n1 and n2 timing
    params.temp_n1_idx = params.n1_idx + params.stim_idx - 1;
    params.temp_loose_n1_idx = params.loose_n1_idx + params.stim_idx - 1;
    params.temp_n2_idx = params.n2_idx + params.stim_idx - 1;
    params.temp_stim_idx = params.stim_indices_0 + params.stim_idx - 1;
    params.temp_tight_stim = params.tight_stim_indices + params.stim_idx-1;
    