% make plots for figure 1 - individual ccep trials separate from average,
% then a representative ccep average that I can manipulate

clear all; close all; clc
addpath(genpath('.'));
%%
%FIGURE_DISPLAY('off');
%% define patient and save directories

locations = cceps_files;
datadir = fullfile(locations.results_folder,'trial_metrics'); 
savedir = fullfile(locations.results_folder,'pub_figs','fig1'); mkdir(savedir);
coords = load('../data/elecs.mat');

%% manually provide indices of representative CCEPs to plot in figure 1
% just plotting both N1 and N2 for simplicity.

plot_guide = struct();

plot_guide.('HUP212').stim = [65];
plot_guide.('HUP212').rec = [69];

for pt = {'HUP212'}
    pt = char(pt);  
    % HUP212: there are some Nans in distance because of missing
    % coordinates for electrode 87 and 88
    %% load data
    % trial level data
    load(fullfile(locations.results_folder,'all_trials',[pt,'_CCEPTrialWaveforms.mat']),'trials');
    load(fullfile(datadir,[pt,'_CCEPTrialSummaryMetrics.mat']));
    
    % load out struct and other patient data, add coordinates
    load([locations.results_folder,'/results_',pt,'_CCEP.mat']);
    clinical = pull_clinical_info([pt,'_CCEP']);    
    
    out = ADD_COORDINATES_BIPOLAR(out,coords);
    out.chLabels_ana = anatomic_location(out.chLabels,clinical,1);
    
    if isempty(out.chLabels_ana)
       out.chLabels_ana = out.chLabels;
    end
    
    % get stim timing
    params = get_stim_timing(out);
    
    %% 
    cmaps = load([locations.data_folder,'colors/mpl_cmaps.mat']);
    max_nt = max(cell2mat(cellfun(@(x) size(x,2),trials(:),'UniformOutput',false))); % maximum number of trials
    %cmap = cmaps.custom_ejc2(round(linspace(1,length(cmaps.custom_ejc2),max_nt),:);
    cmap = spring(max_nt);
    
    waveforms = {'N1','N2'};
    
    % this can be removed after you run compute_trial_summary_metrics.m again
    results.('N1').time = params.n1_time;
    results.('N2').time = params.n2_time;
    
    for wave = waveforms
        wave = char(wave);
        D = results.(wave).D;
        A = results.(wave).A;
        
        good_cceps = ~isempty_c(results.(wave).data) & out.network(1).A>0 & out.network(2).A>0;
        
        %% to what extent do responses change monotonically with increasing trials
        
        [spear_trials,spear_trials_p,spear_trials_p_adj] = COMPUTE_SEQ_SPEARMAN(out,results.(wave).data);
               
        %%
        
        n_plot = min(10,sum(good_cceps,'all'));

        % remove entries where there is no meaningful ccep        
        stim_idx = repmat(1:size(spear_trials,2),size(spear_trials,1),1);
        response_idx = repmat((1:size(spear_trials,1))',1,size(spear_trials,2));
        A = A(good_cceps);
        D = D(good_cceps);
        spear_trials = spear_trials(good_cceps);        
        spear_trials_p_adj = spear_trials_p_adj(good_cceps);
        stim_idx = stim_idx(good_cceps);
        response_idx = response_idx(good_cceps);

        % identify strongest N1s, form table, then sort table
        [~,I] = sort(abs(spear_trials),'descend');                
        info = array2table([stim_idx response_idx spear_trials spear_trials_p_adj INVPRCTILE(A,A) INVPRCTILE(D,D)],'VariableNames',{'Stim','Record','Spear','Spear_p_adj',wave,'Distance'});
        info.StimLabel = out.chLabels_ana(info.Stim);
        info.RecordLabel = out.chLabels_ana(info.Record);    
        
        info = info(I,:);
                            
        for j = 1:length(plot_guide.(pt).stim)
            
            %% plot CCEP waveform
            
            lncol = [40 46 69]/255;
            
            f=figure;
            f=FIGURE_SIZE_CM(f,5,3);        
            tiledlayout(1,1,'padding','compact');
            s_ind = plot_guide.(pt).stim(j);
            r_ind = plot_guide.(pt).rec(j);
            tbl_idx = find(info.Stim == s_ind & info.Record == r_ind);
            
            nt = size(trials{r_ind,s_ind},2);

            % extract EEG data for stim train (stim-recording pair)
            dat_all = trials{r_ind,s_ind}(params.trial_plot_idx,:);
            
            % demean each trial
            %dat = (dat - nanmean(dat))*2;
            %dat_all = (dat_all - nanmean(dat_all))*2;
            
            % get times
            %tw = out.elecs(s_ind).times*1000;
            tw = params.trial_plot_times*1000;
            times = linspace(tw(1),tw(2),size(dat_all,1));
            
            % plot average
            ave_trace = nanmean(dat_all,2);
            %ave_trace() = 1.5*ave_trace(params.temp_n2_idx(1):params.temp_n2_idx(2));
            
            ave=plot(times,ave_trace,'k','LineWidth',1,'Color',lncol);                        
            ave_trace(floor(215e-3*params.fs):floor(380e-3*params.fs)) = ave_trace(floor(215e-3*params.fs):floor(380e-3*params.fs))*1.33;
            hold on;
            ave2 = plot(times+max(times)+400,ave_trace,'k','LineWidth',1,'Color',lncol);                        
            ave_trace(floor(215e-3*params.fs):floor(380e-3*params.fs)) = ave_trace(floor(215e-3*params.fs):floor(380e-3*params.fs))*1.33;
            hold on;
            ave3 = plot(times+2*max(times)+800,ave_trace,'k','LineWidth',1,'Color',lncol);                        
            
            % plot trials just to get y limits
            
            % labels
            ylabel('uV');
            xlabel('Time (ms)');                        
            %xl = xlim; yl= ylim;
            %ylim([yl(1) yl(2)*1.9])
            %xlim([xl(1),xl(2)*6]);
            xticks([0 800 1600]);
            xticklabels([0 1000 2000]);
            prettifyEJC;             
            set(gca,'FontSize',4)         
            
            saveas(f,fullfile(savedir,sprintf('%s_CCEPAverageFig1_%s_S%d_R%d.pdf',wave,pt,s_ind,r_ind)))
                       
            %% plot sin wave
            
            sincol = [62 157 138]/255;
            f=figure;
            f=FIGURE_SIZE_CM(f,5,3);        
            tiledlayout(1,1,'padding','compact');    
            t = linspace(-500,2000,3000);
            plot(t,sin(0.0035*(t+500)),'k','LineWidth',1,'Color',sincol);
            ylim([-1.5 10]);
            xlim([-500 2000]);
            xticks([0 800 1600]);
            xticklabels([0 1000 2000]);
            prettifyEJC;             
            set(gca,'FontSize',4)         
            saveas(f,fullfile(savedir,sprintf('%sSinWave%s_S%d_R%d.pdf',wave,pt,s_ind,r_ind)))
            
            %% plot average and perform piecewise manipulation
                        
        end

        
    end

end







