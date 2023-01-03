function TS_LINE_PLOT_SIMPLE(out,X,sch,chs)

% INPUTS:
% out: out structure
% X: T by N time series
% sch: stim channel index
% chs: recording channel indices
%
% OUTPUTS:
% f: figure handle to plot of e phys time series

nchs = size(X,2);
time = [1:size(X,1)]'/out.other.stim.fs;
sf = 10;
offset_matrix = sf*repmat(1:size(X,2),[size(X,1) 1]);
Xmean = repmat(nanmean(X),[size(X,1) 1]);
Xstd = repmat(nanstd(X),[size(X,1) 1]);
Xz = [X - Xmean] ./ Xstd;

%Xz(isnan(Xz)) = 0; % if any missing trials make them zeros rather than nans to preserve trial order

plot(repmat(time,[1 nchs]),Xz + offset_matrix);
yticklabels(out.chLabels(chs)); yticks([1:nchs]*10);
xlabel('Time (s)'); ylabel('Channel');
