# this script compares number of significant spearman correlations based on method for calculating N1 and N2

rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'eli/miscfxns/packages.R'))
source(paste0(basedir,'eli/miscfxns/miscfxns.R'))
source(paste0(basedir,'eli/miscfxns/statfxns.R'))
source(paste0(basedir,'eli/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'wave_method_spearman/')
dir.create(savedir,recursive = T)

# load distribution of spearman correlations between trial index and N1/N2/pre stim negative control

data <- readMat(paste0(locations$results_folder,'pub_figs/AllPatientsSpearmanLocsSoz.mat'))
spear.results <- unpack.mat(data,'spear.results')
elec.info <- unpack.mat(data,'elec.info')
ave.cceps <- unpack.mat(data,'ave.cceps')

data.methods <- readMat(paste0(locations$results_folder,'pub_figs/AllPatientsWaveQuantMethodsSpearman.mat'))
spear.method.results <- unpack.mat(data.methods,'spear.results.wave.method')

waveforms <- c('N1','N2','PreStim')
df.all <- list()
for(wave in waveforms){
  data.wave <- unpack.mat(spear.method.results,wave)
  methods <- name(names(data.wave))
  p.adj.all.methods <- lapply(methods, function(method) unpack.mat.struct(unpack.mat(data.wave,method),'spear.trials.p.adj'))
  df.plt <- as.data.frame(lapply(p.adj.all.methods,function(p.adj.method)
    sapply(p.adj.method, function(X) sum(X<0.05,na.rm=T))))
  df.plt$pt <- rownames(df.plt)
  
  df.plt <- collapse.columns(df.plt,cnames=methods,groupby = 'pt')
  
  df.plt$wave <- wave
  df.all[[wave]] <- df.plt
  
}

df.plt <- do.call(rbind,df.all)


p <- ggplot(df.plt) + geom_col(aes(x=group,fill=names,y=values),position=position_dodge()) + facet_wrap(~wave,ncol = 1,scales='free') + theme_minimal() +
   ylab('# of Significant Monotonic Trends') + standard_plot_addon() + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=90,vjust=0.5))
ggsave(plot = p,filename = paste0(savedir,'CompareWaveformMethodSpearmans.pdf'),width = 18,height=10,units= 'cm',useDingbats=FALSE)
