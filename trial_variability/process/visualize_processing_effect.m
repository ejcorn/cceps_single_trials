clearvars -except trials trials_raw
%%

sch = 146;
rch = 11;

data_raw = trials_raw{rch,sch};
data = trials{rch,sch};

data_mean = mean(data);
data_sd = std(data);
sf = 10;

f = figure;
subplot(1,2,1);

plot(data_raw./data_sd + sf*repmat(1:size(data,2),[size(data,1) 1]));
title('raw')

subplot(1,2,2);
plot(data./data_sd  + sf*repmat(1:size(data,2),[size(data,1) 1]));
title('processed')

%%
good_cceps = ~isempty_c(trials);

f = @(x) reshape(x,[size(x,1)*size(x,2) 1]);
trials_all = cellfun(f,trials,'UniformOutput',false);
trials_all = [trials(good_cceps)];
trials_all = vertcat(trials_all{:});

f = @(x) reshape(x,[size(x,1)*size(x,2) 1]);
trials_raw_all = cellfun(f,trials_raw,'UniformOutput',false);
trials_raw_all = [trials(good_cceps)];
trials_raw_all = vertcat(trials_raw_all);

%%



