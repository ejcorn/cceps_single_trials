# - for each patient, get # significant rho from recording in each anatomical parcel, stim from every electrode
# - for each patient, get # significant rho from stim in each anatomical parcel, rec from every electrode
# plot across patients and see if there is anatomical consistency in direction of effects

rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'eli/miscfxns/packages.R'))
source(paste0(basedir,'eli/miscfxns/miscfxns.R'))
source(paste0(basedir,'eli/miscfxns/statfxns.R'))
source(paste0(basedir,'eli/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/figS1/')
dir.create(savedir,recursive = T)


# load distribution of spearman correlations between trial index and N1/N2/pre stim negative control

data <- readMat(paste0(locations$results_folder,'pub_figs/AllPatientsSpearmanLocsSoz.mat'))
spear.results <- unpack.mat(data,'spear.results')
elec.info <- unpack.mat(data,'elec.info')
ave.cceps <- unpack.mat(data,'ave.cceps')
brainnetome.names <- read.table(paste0(locations$data_folder,'/nifti/BrainnetomeNodeNames.txt'),header = F,stringsAsFactors = F)$V1
brainnetome.lobes <- read.csv(paste0(locations$data_folder,'/nifti/BNA_subregions_unmerge.csv'),header = T,stringsAsFactors = F)
rownames(brainnetome.lobes) <- brainnetome.names[brainnetome.lobes$Label.ID]
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
  
  # get mean stim and recording rho for each patient for each brainnetome parcel
  bn.rho.stim.pos <- bn.rho.stim.neg <- as.data.frame(matrix(NA,nrow=nparc,ncol=length(pts),dimnames = list(brainnetome.names,pts)))
  bn.rho.rec.pos <- bn.rho.rec.neg <- as.data.frame(matrix(NA,nrow=nparc,ncol=length(pts),dimnames = list(brainnetome.names,pts)))
  
  for(pt in pts){
    df.pt <- df[df$pt == pt,]
    for(roi in 1:nparc){
      
      # for all CCEPs obtained by stimulating a given brainnetome parcel
      # calculate the mean rho
      # test if that mean rho differs from 0
      
      df.pt.roi.stim <- df.pt[df.pt$BN.stim == roi & !is.na(df.pt$BN.stim),]
      if(nrow(df.pt.roi.stim) > 0){
        bn.rho.stim.neg[roi,pt] <- sum(df.pt.roi.stim$p.adj[df.pt.roi.stim$rho<0]<0.05)
        bn.rho.stim.pos[roi,pt] <- sum(df.pt.roi.stim$p.adj[df.pt.roi.stim$rho>0]<0.05)
      }
      
      # for all CCEPs obtained by recording from a given brainnetome parcel
      # calculate the mean rho
      # test if that mean rho differs from 0
      
      df.pt.roi.rec <- df.pt[df.pt$BN.rec == roi & !is.na(df.pt$BN.rec),]
      if(nrow(df.pt.roi.rec) > 0){
        bn.rho.rec.neg[roi,pt] <- sum(df.pt.roi.rec$p.adj[df.pt.roi.rec$rho<0]<0.05)
        bn.rho.rec.pos[roi,pt] <- sum(df.pt.roi.rec$p.adj[df.pt.roi.rec$rho>0]<0.05)
      }
    }
  }
  
  results[[wave]]$Stim$Positive <- bn.rho.stim.pos
  results[[wave]]$Stim$Negative <- bn.rho.stim.neg
  results[[wave]]$Recording$Positive <- bn.rho.rec.pos
  results[[wave]]$Recording$Negative <- bn.rho.rec.neg
  
  df.all[[wave]] <- df
  
}


# plot
for(wave in waveforms){
  #bn.rho.stim.pos * (bn.rho.stim.neg<0.05)
  for(stim.rec in c('Stim','Recording')){
    
    # first plot the number of positive and negative in each region
    bn.rho.pos <- results[[wave]][[stim.rec]]$Positive
    bn.rho.neg <- results[[wave]][[stim.rec]]$Negative
    bn.rho.pos[bn.rho.pos==0 & !is.na(bn.rho.pos)] <- NA
    bn.rho.neg[bn.rho.neg==0 & !is.na(bn.rho.neg)] <- NA
    
    bn.rho.neg <- bn.rho.neg*-1
    
    # only look at regions with effects for at least 2 patients
    multiple.pts <- c(rowSums(!is.na(bn.rho.pos)) + rowSums(!is.na(bn.rho.neg))) >1
    bn.rho.pos <- bn.rho.pos[multiple.pts,]
    bn.rho.neg <- bn.rho.neg[multiple.pts,]
    
    bn.rho.pos$roi.names <- row.names(bn.rho.pos)
    bn.rho.neg$roi.names <- row.names(bn.rho.neg)
    
    df.plt <- collapse.columns(rbind(bn.rho.pos,bn.rho.neg),cnames = pts,groupby = 'roi.names')
    df.plt <- df.plt[!is.na(df.plt$values),]
    df.plt$group <- paste0(df.plt$group,'\n(',brainnetome.lobes[df.plt$group,'Gyrus'],')')
    ylim <- max(abs(df.plt$values))
    
    p <- ggplot(df.plt) + geom_point(aes(x=group,y=values),alpha=0.5,size=1,stroke=0,color=wes_palettes$BottleRocket2[2]) + ggtitle(paste0(wave,': ',stim.rec))+
      scale_y_continuous(limits=c(-ylim-3,ylim+3))+
      theme_bw() + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.title.x = element_blank()) + ylab('# of Significant Effects')+
    standard_plot_addon()

    ggsave(plot = p,filename = paste0(savedir,'CountSignificantRhoPositiveNegative_',stim.rec,wave,'.pdf'),width = 4,height=5.5,units= 'cm',useDingbats=FALSE)  
    
    
  }
}

