rm(list=ls())

basedir <- '/Users/Eli/Dropbox/CNTProjects/CCEPs_Projects/cceps_single_trials/'
source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source('trial_variability/plottingfxns/plottingfxns.R')

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/CCEPMatrices/')
dir.create(savedir,recursive = T)

#####################################################################
### Visualize the CCEP matrix separated by SOZ/Non-SOZ electrodes ###
#####################################################################

data <- readMat(paste0(locations$results_folder,'pub_figs/AllPatientsSpearmanLocsSoz.mat'))
ave.cceps <- unpack.mat(data,'ave.cceps')
spear.results <- unpack.mat(data,'spear.results')

ave.cceps.N1 <- unpack.mat.struct(ave.cceps,'N1')
ave.cceps.N2 <- unpack.mat.struct(ave.cceps,'N2')
pts <- rownames(pt.tbl)

elec.info <- unpack.mat(data,'elec.info')
stim.chs <- unpack.mat.struct(elec.info,'stim.chs')
chLabels.all <- lapply(unpack.mat.struct(elec.info,'chLabels'),cell.to.vec)
rec.chs <- unpack.mat.struct(elec.info,'rec.chs')
fs <- unpack.mat.struct(elec.info,'fs')

# tile soz indicator to matrix size so you know if SOZ was stimulated for each CCEP
soz.all.pts <- unpack.mat.struct(elec.info,'soz')
#soz.all.pts <- lapply(soz.all.pts, function(x) ifelse(x,yes='SOZ',no='Non-SOZ'))
#soz.all.pts.rec <- lapply(soz.all.pts, function(x) matrix(x,nrow=length(x),ncol=length(x))) # check: colSums(soz.all.pts.rec$HUP211=='SOZ') should be 25 for every electrode
#soz.all.pts.stim <- lapply(soz.all.pts, function(x) t(matrix(x,nrow=length(x),ncol=length(x)))) # check: rowSums(soz.all.pts.stim$HUP211=='SOZ') should be 25 for every electrode

pts <- names(soz.all.pts)
pts[!is.na(sapply(soz.all.pts,sum))]

p.list <- df.lines <- list()
for(pt in pts){
  soz <- soz.all.pts[[pt]]
  #imagesc(good.cceps[order(soz),order(soz)])
  
  X <- ave.cceps.N1[[pt]][order(soz),order(soz)]
  rownames(X) <- chLabels.all[[pt]][order(soz)]
  colnames(X) <- chLabels.all[[pt]][order(soz)]
  
  no.stim.mask <- colMeans(is.na(X))<1
  X <- X[,no.stim.mask]
  
  
  p <- imagesc(X,clim=c(0,40)) + coord_equal() + x_text_90() + theme(axis.text = element_text(size=4))+
    ggtitle(paste(chLabels.all[[pt]][which(soz==1)],collapse=', ')) 
  
  ggsave(plot = p,filename = paste0(savedir,pt,'N1Matrix.pdf'),width = 24,height=24,units= 'cm',useDingbats=FALSE)
  
  # make plot with lines over SOZ and no tick labels to make a nice figure
  soz.ord <- soz[order(soz)]
  soz.line.x <- 2.5+length(soz)-min(which(soz.ord==1))
  soz.line.y <- -1.5+min(which(soz.ord[no.stim.mask]==1))
  if(soz.line.y==Inf){soz.line.y <- NA}
  df.lines[[pt]] <- data.frame(x=soz.line.x,y=soz.line.y)
  p.list[[pt]] <- imagesc(X,clim=c(0,30)) +
    geom_hline(data = df.lines[[pt]],aes(yintercept=x),color='white',size=0.25)+
    geom_vline(data=df.lines[[pt]],aes(xintercept=y),color='white',size=0.25)+
    theme(axis.text = element_blank(),axis.ticks = element_blank())+
    ggtitle(pt) + nice_cbar()
  
}

pts.plot <- c('HUP213','HUP214','HUP218','HUP219','HUP223','HUP225')
p.all <- plot_grid(plotlist = p.list[pts.plot],align = 'hv',nrow=2)
ggsave(plot = p.all,filename = paste0(savedir,'AllPatientsN1MatricesBySOZ.pdf'),width = 18,height=12,units= 'cm',useDingbats=FALSE)
 