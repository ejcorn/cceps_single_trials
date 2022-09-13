clearvars -except trials trials_raw
%%

good_cceps = ~isempty_c(trials);

% function to concatenate trials together
f = @(x) reshape(x,[size(x,1)*size(x,2) 1]); 
trials_all = cellfun(f,trials,'UniformOutput',false);
[~,rchs] = find(good_cceps);
for rch = unique(rchs)'
    trials_rch = trials_all(rch,:);
end

% I need to go back and download time series data for all sch-rch pairs