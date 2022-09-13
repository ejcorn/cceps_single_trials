rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig4/')
dir.create(savedir,recursive = T)


fnames <- Sys.glob(paste0(locations$results_folder,'/powercalc_peakfreq_AL/peaks_*.csv'))
files <- lapply(fnames,function(x) read.csv(x,header=F))
x <- do.call(rbind,files)
x$retained <- ifelse(x$V1 <14,yes='Retained',no='Discarded')

p <- ggplot() + geom_histogram(aes(x=x$V1,fill=x$retained),alpha=0.5) + geom_vline(xintercept = 14,linetype='dashed') + theme_classic() +
  scale_y_continuous(expand=c(0,0)) + 
  scale_fill_manual(values=c(wes_palettes$Rushmore[3],wes_palettes$BottleRocket1[3]),breaks=c('Retained','Discarded')) +
  xlab('Lowest Frequency Peak (Hz)') + ylab('# of Electrodes') + ggtitle('Peak Selection') + standard_plot_addon() +
  theme(legend.title = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key.size = unit(0.1,'cm')) 
p
ggsave(plot = p,filename = paste0(savedir,'Fig4a_bottom.pdf'),width = 4,height=4,units= 'cm')
