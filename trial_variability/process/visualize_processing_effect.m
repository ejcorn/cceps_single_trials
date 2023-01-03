%clearvars -except trials trials_raw
%%
pt = 'HUP219';
good_cceps = out.network(1).A> 10 & out.network(2).A> 10;
[row,col] = find(good_cceps);
j = 106;
sch = col(j);
rch = row(j);

data_raw = trials_raw{rch,sch};
data = trials{rch,sch};

data_mean = mean(data);
data_sd = std(data);
sf = 10;

f = figure;
subplot(1,3,1);

plot(data_raw./data_sd + sf*repmat(1:size(data,2),[size(data,1) 1]));
xlim([300 900])
title('raw')

subplot(1,3,2);
plot(data./data_sd  + sf*repmat(1:size(data,2),[size(data,1) 1]));
xlim([300 900])
title('processed')

data_alt = PREPROCESS_CCEPS_LINE_ONLY(out,data_raw',true)';
subplot(1,3,3);
plot(data_alt./data_sd  + sf*repmat(1:size(data,2),[size(data,1) 1]));
xlim([300 900])
title('processed - line noise filter only')
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



