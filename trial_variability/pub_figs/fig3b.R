# - for each patient, get # significant rho from recording in each anatomical parcel, stim from every electrode
# - for each patient, get # significant rho from stim in each anatomical parcel, rec from every electrode
# plot within each patient the nmber of times positive and negative effects overlap

rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/miscfxns/statfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig3/fig3b/')
dir.create(savedir,recursive = T)

# load distribution of spearman correlations between trial index and N1/N2/pre stim negative control

data <- readMat(paste0(locations$results_folder,'pub_figs/AllPatientsSpearmanLocsSoz.mat'))
spear.results <- unpack.mat(data,'spear.results')
elec.info <- unpack.mat(data,'elec.info')
ave.cceps <- unpack.mat(data,'ave.cceps')
brainnetome.names <- read.table(paste0(locations$data_folder,'/nifti/BrainnetomeNodeNames.txt'),header = F,stringsAsFactors = F)$V1
nparc <- length(brainnetome.names)

# loop through wave forms, look at number of significant spearman effects for SOZ and Non-SOZ electrodes

waveforms <- c('N1','N2')
df.all <- list()
results <- list()
wave <- 'N1'
pt <- 'HUP223'

for(wave in waveforms){
  
  # unpack mat struct
  data.wave <- unpack.mat(spear.results,wave)
  spear.all.pts <- unpack.mat.struct(data.wave,'spear.trials')
  good.cceps.all.pts <- unpack.mat.struct(data.wave,'good.cceps')
  p.adj.all.pts <- unpack.mat.struct(data.wave,'spear.trials.p.adj')
  
  D.all.pts <- unpack.mat.struct(elec.info,'D')
  chLabels.all <- lapply(unpack.mat.struct(elec.info,'chLabels'),cell.to.vec)
  chLabels.ana.all <- lapply(unpack.mat.struct(elec.info,'chLabels.ana'),cell.to.vec)
  stim.order.all <- lapply(unpack.mat.struct(elec.info,'chLabel.stim.order'),cell.to.vec)
  
  Brainnetome.all.pts <- unpack.mat.struct(elec.info,'Brainnetome')
  
  # matrix whose elements correspond to Brainnetome parcel assignment for stim or recording electrode of each CCEP
  Brainnetome.all.pts.stim <- lapply(Brainnetome.all.pts, function(x) t(matrix(x,nrow=length(x),ncol=length(x)))) # check: colSums(soz.all.pts.stim$HUP211=='SOZ') should be 25 for every electrode
  Brainnetome.all.pts.rec <- lapply(Brainnetome.all.pts, function(x) matrix(x,nrow=length(x),ncol=length(x))) # check: rowSums(soz.all.pts.rec$HUP211=='SOZ') should be 25 for every electrode
  
  # matrix whose elements correspond to name of stimulating electrode for each matrix element (i.e. CCEP)
  rec.name <- lapply(chLabels.all, function(x) matrix(x,nrow=length(x),ncol=length(x)))
  rec.name.ana <- lapply(chLabels.ana.all, function(x) matrix(x,nrow=length(x),ncol=length(x)))
  stim.name <- lapply(chLabels.all, function(x) t(matrix(x,nrow=length(x),ncol=length(x))))
  stim.name.ana <- lapply(chLabels.ana.all, function(x) t(matrix(x,nrow=length(x),ncol=length(x))))
  
  pts <- name(names(spear.all.pts))
  pt.indicator <- lapply(pts, function(pt) matrix(pt,nrow=nrow(spear.all.pts[[pt]]),ncol=ncol(spear.all.pts[[pt]])))
  
  df <- data.frame(rho=list.mat.to.vec(spear.all.pts),
                   g=list.mat.to.vec(good.cceps.all.pts),
                   p.adj=list.mat.to.vec(p.adj.all.pts),
                   BN.stim=list.mat.to.vec(Brainnetome.all.pts.stim),
                   BN.rec=list.mat.to.vec(Brainnetome.all.pts.rec),
                   D=list.mat.to.vec(D.all.pts),
                   stim.name=list.mat.to.vec(stim.name),
                   stim.name.ana=list.mat.to.vec(stim.name.ana),
                   rec.name=list.mat.to.vec(rec.name),
                   rec.name.ana=list.mat.to.vec(rec.name.ana),
                   pt=list.mat.to.vec(pt.indicator),
                   stringsAsFactors = F)
  df <- df[df$g==1,]
  
  
  for(pt in pts){
    
    chLabels.pt <- chLabels.all[[pt]]
    nchs.pt <- length(chLabels.pt)
    
    elec.sig.stim.pos <- elec.sig.stim.neg <- as.data.frame(matrix(NA,nrow=nchs.pt,ncol=1,dimnames = list(chLabels.pt,NULL)))
    elec.sig.rec.pos <- elec.sig.rec.neg <- as.data.frame(matrix(NA,nrow=nchs.pt,ncol=1,dimnames = list(chLabels.pt,NULL)))
    df.pt <- df[df$pt == pt,]
    for(roi in chLabels.pt){
      
      # for all CCEPs obtained by stimulating a given brainnetome parcel
      # calculate the mean rho
      # test if that mean rho differs from 0
      
      df.pt.elec.stim <- df.pt[df.pt$stim.name == roi,]
      if(nrow(df.pt.elec.stim) > 0){
        elec.sig.stim.neg[roi,] <- sum(df.pt.elec.stim$p.adj[df.pt.elec.stim$rho<0]<0.05)
        elec.sig.stim.pos[roi,] <- sum(df.pt.elec.stim$p.adj[df.pt.elec.stim$rho>0]<0.05)
      }
      
      # for all CCEPs obtained by recording from a given brainnetome parcel
      # calculate the mean rho
      # test if that mean rho differs from 0
      
      df.pt.elec.rec <- df.pt[df.pt$rec.name == roi,]
      if(nrow(df.pt.elec.rec) > 0){
        elec.sig.rec.neg[roi,] <- sum(df.pt.elec.rec$p.adj[df.pt.elec.rec$rho<0]<0.05)
        elec.sig.rec.pos[roi,] <- sum(df.pt.elec.rec$p.adj[df.pt.elec.rec$rho>0]<0.05)
      }
      
    }
    
    results[[wave]]$Stim$Positive[[pt]] <- elec.sig.stim.pos
    results[[wave]]$Stim$Negative[[pt]] <- elec.sig.stim.neg
    results[[wave]]$Recording$Positive[[pt]] <- elec.sig.rec.pos
    results[[wave]]$Recording$Negative[[pt]] <- elec.sig.rec.neg
    
  }
  
  
  
  df.all[[wave]] <- df
  
}


# plot
for(wave in waveforms){
  #elec.sig.stim.pos * (elec.sig.stim.neg<0.05)
  for(stim.rec in c('Stim','Recording')){
    
    
    elec.sig.pos <- results[[wave]][[stim.rec]]$Positive
    elec.sig.neg <- results[[wave]][[stim.rec]]$Negative
    
    df.plt <- list()
    for(pt in pts){
      elec.sig.pos[[pt]][elec.sig.pos[[pt]]==0 & !is.na(elec.sig.pos[[pt]])] <- NA
      elec.sig.neg[[pt]][elec.sig.neg[[pt]]==0 & !is.na(elec.sig.neg[[pt]])] <- NA
      
      elec.sig.and <- elec.sig.neg[[pt]] & elec.sig.pos[[pt]] # get ROIs with both positive and negative effects in same patient
      elec.sig.neg[[pt]][elec.sig.and] <- NA # remove electrodes with both positive and negative effects
      elec.sig.pos[[pt]][elec.sig.and] <- NA # remove electrodes with both positive and negative effects
      
      elec.sig.neg.pt <- data.frame(roi.count=sum(elec.sig.neg[[pt]],na.rm=T),pt=pt,and.or='Negative')
      elec.sig.pos.pt <- data.frame(roi.count=sum(elec.sig.pos[[pt]],na.rm=T),pt=pt,and.or='Positive')
      elec.sig.and <- data.frame(roi.count=sum(elec.sig.and,na.rm=T),pt=pt,and.or='Bidirectional')
      df.plt[[pt]] <- rbind(elec.sig.neg.pt,elec.sig.pos.pt,elec.sig.and)
    }
    
    df.plt <- do.call(rbind,df.plt)
    
    p <- ggplot(df.plt) + geom_col(aes(x=pt,y=roi.count,fill=and.or),position=position_dodge()) + theme_bw() +
      ylab('# of Electrodes with\nSignificant Effects') + 
      scale_fill_manual(values=c(wes_palettes$GrandBudapest2[4],wes_palettes$BottleRocket2[2:3]))+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),axis.title.x = element_blank()) + ggtitle(paste0(wave,': ',stim.rec)) +
      standard_plot_addon() + theme(legend.title = element_blank(),legend.key.size = unit(0.1,'cm'))
    ggsave(plot = p,filename = paste0(savedir,'CountUnidirectionalBidirectionalElectrodes_',stim.rec,wave,'.pdf'),width = 6,height=4,units= 'cm',useDingbats=FALSE)  
    
    
  }
}

