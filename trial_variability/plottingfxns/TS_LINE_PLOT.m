function TS_LINE_PLOT(out,X,sch,k,sch_plot,plot_all)

% INPUTS:
% out: out structure
% X: T by N time series
% sch: stim channel
% k: number of electrodes to plot (default 3)
% sch_plot: plot k electrodes centered around sch_plot (default sch)
%
% OUTPUTS:
% f: figure handle to plot of e phys time series
%

if ~exist('k','var')
    k=3;
end

if ~exist('sch_plot','var')
    sch_plot = sch;
end

if ~exist('plot_all','var')
    plot_all = '';
end

%{
subplot(1,3,1);
plot(zscore(X) + [1:size(X,2)])
subplot(1,3,2);
plot(X(:,sch));
subplot(1,3,3);
%}
if ~strcmp(plot_all,'plot all')
    chs = [(sch_plot-k):(sch_plot+k)];
    nchs = (2*k+1);
else
    nchs = size(X,2);
    chs = 1:nchs;    
end

time = [1:size(X,1)]'/out.other.stim.fs;
plot(repmat(time,[1 nchs]),zscore(X(:,chs)) + [1:nchs]*10);
yticklabels(out.chLabels(chs)); yticks([1:nchs]*10); title(sprintf('Stim Electrode %s',out.chLabels{sch}));
xlabel('Time (s)'); ylabel('Channel');
