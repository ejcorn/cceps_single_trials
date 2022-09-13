function p = CORR_SPEAR_P(x,y)

% compute spearman correlation and return p-value

[~,p] = corr(x,y,'rows','pairwise','type','spearman');