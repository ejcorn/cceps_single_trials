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

pt = 'HUP212';
load([locations.results_folder,'/results_',pt,'_CCEP.mat']);
clinical = pull_clinical_info([pt,'_CCEP']);    
out = ADD_COORDINATES_BIPOLAR(out,coords);
out.chLabels_ana = anatomic_location(out.chLabels,clinical,1);

if isempty(out.chLabels_ana)
   out.chLabels_ana = out.chLabels;
end

stim_chs = find(out.stim_chs);
[~,data] = get_stim_trials_sch(out,stim_chs(1));

%%
good_cceps = out.network(1).A>0 & out.network(2).A>0;
X = data.values(:,40:44);
load(fullfile(locations.data_folder,'colors/mpl_cmaps.mat'))

f=figure;
f=FIGURE_SIZE_CM(f,8,2);
tiledlayout(1,1,'padding','compact');    
plot(zscore(X) + 10*[1:size(X,2)])
xticklabels([]); yticklabels([]);

set(gca,'LineWidth',1);
set(gca,'TickLength',[0 0]);
set(gca,'FontName','arial')
xlim([0 length(X)])

cmap = Blues;
cmap = cmap(floor(linspace(40,length(cmap),5)),:);
colororder(cmap)
set(gca,'FontSize',4)         
saveas(f,fullfile(locations.results_folder,'pub_figs','fig1','TSExample.pdf'));
