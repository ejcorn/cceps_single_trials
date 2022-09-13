% To study rostral and caudal hippocampus, look at just one electrode in
% each ROI
% ideally pick the electrode with the greatest theta power peak
% for now i will just test all the electrodes separately

% using the circ_corrcl function from circstat function 


%%
clear all; close all; clc
addpath(genpath('.'));
%%
FIGURE_DISPLAY('off');
%% define patient and save directories

locations = cceps_files;
datadir = fullfile(locations.results_folder,'trial_metrics'); 
savedir = fullfile(locations.results_folder,'hippocampal_theta_dependence_circ_corrcl'); mkdir(savedir);
savedir_R = fullfile(locations.results_folder,'pub_figs','fig4'); mkdir(savedir_R);
coords = load('../data/elecs.mat');

%% Set other parameters

% indices of rostral and caudal hippocampus in the Brainnetome atlas
ROIs = 215:218;

% define waveform names
waveforms = {'N1','N2','PreStim'};

% out struct
clcorr_results = struct();

%%
for pt = locations.subjects
    
    pt = char(pt);  
    % HUP212: there are some Nans in distance because of missing
    % coordinates for electrode 87 and 88
    %% load data
    
    % CCEP trial EEG data
    load(fullfile(locations.results_folder,'all_trials',[pt,'_CCEPTrialWaveforms.mat']),'trials');
    
    % wave form quantification    
    load(fullfile(locations.results_folder,'trial_metrics',[pt,'_CCEPTrialSummaryMetrics.mat']));
    
    % power/phase data
    %features = load(fullfile(locations.results_folder,'powercalc_peakfreq_AL',sprintf('%s_CCEPPrestimFeatures_PeakFreq.mat',pt)));
    features = load(fullfile(locations.results_folder,'prestim_peakfreq_EJC/',sprintf('%s_CCEPPrestimFeatures_PeakFreq_0.005.mat',pt)));
    fooof_peaks = csvread(fullfile(locations.results_folder,'powercalc_peakfreq_AL',sprintf('peaks_%s.csv',pt)));
    
    % load out struct and other patient data, add coordinates
    load([locations.results_folder,'/results_',pt,'_CCEP.mat']);
    clinical = pull_clinical_info([pt,'_CCEP']);    
    
    out = ADD_COORDINATES_BIPOLAR(out,coords);
    out.chLabels_ana = anatomic_location(out.chLabels,clinical,1);
    
    if isempty(out.chLabels_ana)
        %
        out.chLabels_ana = coords.info(find(strcmp({coords.info.name},pt))).elecs.anatomy;
        out.chLabels_ana(isempty_c(out.chLabels_ana)) = {''};

        % keeping this in here
        out.chLabels_ana = [repmat({'ECG EEG'},[8 1]);out.chLabels_ana];

        if isempty(out.chLabels_ana) | length(out.chLabels) ~= length(out.chLabels_ana)
            out.chLabels_ana = out.chLabels;       
        end
    end
    %}
    % get stim timing
    params = get_stim_timing(out);           
    
  
    feat_elec_indices = find(ismember(out.parcellation.elec_parcels.Brainnetome,ROIs))';    
    feat_elec_indices = feat_elec_indices(fooof_peaks(feat_elec_indices) < 15); % don't analyze any electrodes if they are missing slow oscillations
    
    % Loop through ROIs (anterior and posterior hippocampus
    for feat_elec = feat_elec_indices
        
        % get mask of which electrodes belong to parcels of interest in
        % Brainnetome atlas

            % Loop through features (high theta, low theta) - only do phase
            % for now
            for j_feat = find(strcmp('inst_phase',features.data.feat_name)) %1:length(features.data.feat_name)

                % generate label for ROI and feature combo i.e. rostral
                % hippocampal high theta instantaneous phase
                roi_feat_label = sprintf('%s_%s',features.data.feat_name{j_feat},features.data.band_names{j_feat});                
                struct_roi_feat_label = sprintf('%s_%s',out.chLabels{feat_elec},features.data.feat_name{j_feat});
                
                savedir_j = fullfile(savedir,roi_feat_label);
                mkdir(savedir_j);

                % for the jth feature, select only columns corresponding to
                % that feature
                feat_cell = features.data.features(:,j_feat);

                % Average feature over target ROIs only
                any_stim = ~isempty_c(feat_cell);                         

                % Isolate feature for SINGLE ELECTRODE belonging to the
                % specified ROIs
                feat_cell(any_stim) = cellfun(@(x) x(feat_elec,:)',feat_cell(any_stim),'UniformOutput',false);

                % tile feature matrix so that columns contain the feature of
                % interest averaged over the target ROIs, during stimulation of
                % the electrode on the column index
                % i.e. average anterior hippocampal theta phase while 
                feat_cell = repmat(feat_cell',[length(feat_cell) 1]);

                for wave = waveforms
                    wave = char(wave);
                    D = results.(wave).D;
                    A = results.(wave).A;

                    good_cceps = ~isempty_c(results.(wave).data) & out.network(1).A>0 & out.network(2).A>0;

                    %% to what extent do responses change monotonically with pre stim oscillatory features

                    [clcorr_trials,clcorr_trials_p,clcorr_trials_p_adj] = COMPUTE_TRIAL_CIRC_CORRCL(out,feat_cell,results.(wave).data);

                    clcorr_results.(wave).(pt).(struct_roi_feat_label).clcorr_trials = clcorr_trials;
                    clcorr_results.(wave).(pt).(struct_roi_feat_label).clcorr_trials_p_adj = clcorr_trials_p_adj;
                    %%

                    n_plot = min(10,sum(good_cceps,'all'));

                    % remove entries where there is no meaningful ccep        
                    stim_idx = repmat(1:size(clcorr_trials,2),size(clcorr_trials,1),1);
                    response_idx = repmat((1:size(clcorr_trials,1))',1,size(clcorr_trials,2));
                    A = A(good_cceps);
                    D = D(good_cceps);
                    clcorr_trials = clcorr_trials(good_cceps);        
                    clcorr_trials_p_adj = clcorr_trials_p_adj(good_cceps);
                    stim_idx = stim_idx(good_cceps);
                    response_idx = response_idx(good_cceps);

                    % identify strongest N1s, form table, then sort table
                    [~,I] = sort(abs(clcorr_trials),'descend');                
                    info = array2table([stim_idx response_idx clcorr_trials clcorr_trials_p_adj INVPRCTILE(A,A) INVPRCTILE(D,D)],'VariableNames',{'Stim','Record','CLcorr','CLcorr_p_adj',wave,'Distance'});
                    info.StimLabel = out.chLabels_ana(info.Stim);
                    info.RecordLabel = out.chLabels_ana(info.Record);    
                    info.Patient = repmat(pt,[size(info,1) 1]);
                    
                    info = info(I,:);
                    
                    % exclude stim-recording pairs that can't be
                    % calculated, should be when you are stimulating at the
                    % ROI for theta and there are no other electrodes in
                    % the ROI with good signal during stim
                    
                    info = info(~isnan(info.CLcorr),:);
%{
                    f=figure;
                    [sp1,sp2] = subplot_ind2(n_plot);

                    for j = 1:n_plot
                        s_ind = info.Stim(j);
                        r_ind = info.Record(j);

                        subplot(sp1,sp2,j);

                        %nw = plot(feat_cell{r_ind,s_ind},results.(wave).data{r_ind,s_ind},'.','MarkerSize',5,'Color',[1 0 0]);
                        nw = polarplot(feat_cell{r_ind,s_ind},results.(wave).data{r_ind,s_ind},'.','MarkerSize',5,'Color',[1 0 0]);
                        pax=gca;
                        thetaticks(-180:45:180);
                        thetalim([-180 180]);
                        pax.ThetaAxisUnits = 'radians';

                        % plot peak window (n1, n2) as shaded rectangle 
                        title({info.StimLabel{j},['->',info.RecordLabel{j}],sprintf('%s %0.3g Hz',features.data.feat_name{j_feat},fooof_peaks(feat_elec))}, 'Interpreter', 'none');
                        legend([nw],{sprintf('%s: %0.1fth %%ile\nDistance: %0.1fth %%ile\nRho = %.2f\np= %.1d',wave,info.(wave)(j),info.Distance(j),info.CLcorr(j),info.CLcorr_p_adj(j))},'Location','southoutside');
                        prettifyEJC;                        

                    end

                    f=FIGURE_SIZE_CM(f,24,16);
                    saveas(f,fullfile(savedir_j,[pt,'_PlotTopCLcorr',out.chLabels{feat_elec},'FeatureVs',wave,'.pdf']));
                    close(f);
                    
%{                    
                    f=figure;
                    [sp1,sp2] = subplot_ind2(n_plot);

                    for j = 1:n_plot
                        s_ind = info.Stim(j);
                        r_ind = info.Record(j);

                        subplot(sp1,sp2,j);
                        
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

                        % plot average            
                        ave=plot(times,nanmean(dat_all,2),'k','LineWidth',0.5); hold on;
                        % plot individual trials
                        trls=plot(repmat(times',[1 size(dat,2)]),dat,'LineWidth',1); colororder(cmap); hold on;            
                        % plot peak window (n1, n2) as shaded rectangle           
                        hold on;            
                        yl = get(gca,'YLim');
                        re=rectangle('Position',[1000*results.(wave).time(1) -10000 1000*diff(results.(wave).time) 50000],...
                            'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5 0]);
                        set(gca,'YLim',yl);
                        % re order plot elements
                        set(gca,'Children',[ave ;trls;re]);

                        % labels
                        sl = info.StimLabel{j};
                        rl = info.RecordLabel{j};
                        title({[sl,'->'],rl});
                        ylabel('uV');
                        xlabel('Time (ms)');                        

                        legend([ave],{sprintf('%s: %0.1fth %%ile\nDistance: %0.1fth %%ile\nRho = %.2f\np= %.1d',wave,info.(wave)(j),info.Distance(j),info.CLcorr(j),info.CLcorr_p_adj(j))},'Location','southoutside');
                        
                        prettifyEJC;             
                        set(gca,'FontSize',6)
                        
                        
                    end
                    
                    f=FIGURE_SIZE_CM(f,24,16);
                    saveas(f,fullfile(savedir_j,[pt,'_PlotCCEPTraceTopCLcorr',out.chLabels{feat_elec},'FeatureVs',wave,'.pdf']));
                    close(f);
%}
                    fprintf('%s %s %s %0.3g Hz: %d significant\n',pt,wave,struct_roi_feat_label,fooof_peaks(feat_elec),sum(clcorr_trials_p_adj<0.05))
                    
                end
            end
        end
end

save(fullfile(savedir_R,'AllPatientsHippocampalPhaseCLcorr.mat'),'clcorr_results')

FIGURE_DISPLAY('on');