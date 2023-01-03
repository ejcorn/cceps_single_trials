rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/miscfxns/statfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig5/')
dir.create(savedir,recursive = T)

load(file=paste0(savedir,'fig5cde/fig5cdeSOZOutcomeModels.RData'))
# patient - N1 ccep - N1 std - N1 rho - N2 ccep - N2 std - N2 rho

pts <- names(m.soz.pred$N1$std)

pt.tbl <- data.frame(row.names = pts,check.names = F)

waves <- c('N1','N2')
cfg.all <- list(list(
  l.name = 'ccep',
  c.name = 'scale(log10(ccep))',
  nice.name = 'CCEP Amp. (OR, 95% CI)'
),
list(
  l.name = 'std',
  c.name = 'scale(std.trials)',
  nice.name = 'Std. (OR, 95% CI)'
),
list(
  l.name = 'rho',
  c.name = 'scale(rho)',
  nice.name = 'Rho (OR, 95% CI)'
)
)

for(wave in waves){
  for(cfg in cfg.all){
    m.all <- m.soz.pred[[wave]][[cfg$l.name]]
    #m.all <- lapply(m.all,lm.beta)
    pt.tbl[,paste0(wave,' - ',cfg$nice.name)] <- sapply(m.all,function(m) get.beta.OR.95ci(m,cfg$c.name))
    pvals <- p.adjust(sapply(m.all,function(m) get.coef.pval(m,cfg$c.name)),method='fdr')
    pt.tbl[,paste0(wave,' - ',cfg$nice.name)] <- paste(pt.tbl[,paste0(wave,' - ',cfg$nice.name)],p.signif(pvals,ns = ''))
  }
}

xtable(pt.tbl,caption = 'Odds ratios for CCEP metrics in estimating which CCEPs belong to the SOZ.',label = 'table:tableS2')

