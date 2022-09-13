function PLOT_TRIALS_WAVE_WINDOW(params,eeg,window)

nt = size(eeg,2);
nt_plot = 30; % only show subset of trials to clearly visualize trend
cmap = spring(nt_plot);

% extract EEG data for stim train (stim-recording pair)
dat_all = eeg;
%trials_plot_idx = floor(linspace(1,nt,nt_plot));
dat = dat_all; % only show subset of trials

% demean each trial
%dat = (dat - nanmean(dat))*2;
%dat_all = (dat_all - nanmean(dat_all))*2;

% get times
%tw = out.elecs(s_ind).times*1000;
tw = params.time_to_take*1000;
times = linspace(tw(1),tw(2),size(dat,1));

% plot individual trials
trls=plot(repmat(times',[1 size(dat,2)]),dat,'LineWidth',1); 
colororder(cmap); hold on;            

% plot average            
ave=plot(times,nanmean(dat_all,2),'k','LineWidth',0.5); hold on;

% plot peak window (n1, n2) as shaded rectangle           
hold on;            
yl = get(gca,'YLim');
re=rectangle('Position',[1000*window(1) -10000 1000*diff(window) 50000],...
    'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5 0]);
set(gca,'YLim',yl);
% re order plot elements
set(gca,'Children',[ave ;trls;re]);

% labels
ylabel('uV');
xlabel('Time (ms)');                        

prettifyEJC;             
set(gca,'FontSize',6)