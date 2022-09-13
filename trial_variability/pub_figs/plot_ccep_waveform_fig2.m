% plots for figure 2 - showing average CCEP overlaid with individual
% trials, and waveform quantities scatterplotted against trial index

clear all; close all; clc
addpath(genpath('.'));
%%
FIGURE_DISPLAY('off');
%% define patient and save directories

locations = cceps_files;
datadir = fullfile(locations.results_folder,'trial_metrics'); 
savedir = fullfile(locations.results_folder,'pub_figs','fig2'); mkdir(savedir);
coords = load('../data/elecs.mat');

%% manually provide indices of representative CCEPs to plot in figures
% just plotting both N1 and N2 for simplicity.

plot_guide = struct();

plot_guide.('HUP212').stim = [7];
plot_guide.('HUP212').rec = [19];

plot_guide.('HUP213').stim = [10];
plot_guide.('HUP213').rec = [137];

plot_guide.('HUP222').stim = [167];
plot_guide.('HUP222').rec = [188];

plot_guide.('HUP225').stim = [37];
plot_guide.('HUP225').rec = [45];

for pt = {'HUP213','HUP212','HUP225','HUP222'}
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
            
            f=figure;
            f=FIGURE_SIZE_CM(f,4.25,4);        
            tiledlayout(1,1,'padding','compact');
            s_ind = plot_guide.(pt).stim(j);
            r_ind = plot_guide.(pt).rec(j);
            tbl_idx = find(info.Stim == s_ind & info.Record == r_ind);
            
            nt = size(trials{r_ind,s_ind},2);
            nt_plot = 3; % only show subset of trials to clearly visualize trend
            cmap = spring(nt_plot);

            % extract EEG data for stim train (stim-recording pair)
            dat_all = trials{r_ind,s_ind}(params.trial_plot_idx,:);
            trials_plot_idx = floor(linspace(1,nt,nt_plot));
            dat = dat_all(:,trials_plot_idx); % only show subset of trials
            
            % demean each trial
            %dat = (dat - nanmean(dat))*2;
            %dat_all = (dat_all - nanmean(dat_all))*2;
            
            % get times
            %tw = out.elecs(s_ind).times*1000;
            tw = params.trial_plot_times*1000;
            times = linspace(tw(1),tw(2),size(dat,1));
            
            % plot individual trials
            trls=plot(repmat(times',[1 size(dat,2)]),dat,'LineWidth',1); 
            colororder(cmap); hold on;            
            
            % plot average            
            ave=plot(times,nanmean(dat_all,2),'k','LineWidth',0.5); hold on;
            
            % plot peak window (n1, n2) as shaded rectangle           
            hold on;            
            yl = get(gca,'YLim');
            re=rectangle('Position',[1000*results.(wave).time(1) -10000 1000*diff(results.(wave).time) 50000],...
                'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5 0]);
            set(gca,'YLim',yl);
            % re order plot elements
            set(gca,'Children',[ave ;trls;re]);
            
            % labels
            sl = info.StimLabel{tbl_idx};
            rl = info.RecordLabel{tbl_idx};
            title(['Stim: ',sl,'   Record:',rl]);
            ylabel('uV');
            xlabel('Time (ms)');                        
            
            prettifyEJC;             
            set(gca,'FontSize',6)
            saveas(f,fullfile(savedir,sprintf('%s_CCEPAverageAndIndividualTrialsFig2_%s_S%d_R%d.pdf',wave,pt,s_ind,r_ind)))
            close(f);
            
            %% plot colormap
            
            f=figure;
            f=FIGURE_SIZE_CM(f,4.25,4);
            caxis([1 nt]); colormap('spring'); cb = colorbar('southoutside','Ticks',trials_plot_idx);            
            prettifyEJC;             
            set(gca,'FontSize',6)         
            saveas(f,fullfile(savedir,sprintf('Fig2Colormap.pdf',wave,pt,s_ind,r_ind)))
            close(f);
                        
            %% plot N1/N2 strength as scatter against trials
            
            f=figure;
            f=FIGURE_SIZE_CM(f,4,4);
            nw = scatter(1:nt,results.(wave).data{r_ind,s_ind},10,'filled',...
                'MarkerEdgeColor','none',...[1 0.5 0.3],...
                'MarkerFaceColor',[1 0.5 0.3],...
                'Markerfacealpha',0.5);            
            text(0.02,0.93,['\rho',sprintf(' = %.2f\np_{FDR}= ',info.Spear(tbl_idx)),SCINOT_P(info.Spear_p_adj(tbl_idx))],...
                'Units','normalized','FontSize',4);
            xlabel('Trial Number');
            ylabel([sl,' -> ',rl,' ',wave, ' (uV)']);
            prettifyEJC;
            set(gca,'FontSize',6)
            saveas(f,fullfile(savedir,sprintf('%s_Vs_TrialFig2_%s_S%d_R%d.pdf',wave,pt,s_ind,r_ind)))
            close(f);
                        
        end

        
    end

end
FIGURE_DISPLAY('on');