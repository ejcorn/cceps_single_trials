# test for relationship between CCEP amplitude, Rho, Fail rate and amount of spikes in each stimulating or recording
# electrode, using linear model with permutation test to shuffle the electrodes and generate null distribution of betas

rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'eli/miscfxns/packages.R'))
source(paste0(basedir,'eli/miscfxns/miscfxns.R'))
source(paste0(basedir,'eli/miscfxns/statfxns.R'))
source(paste0(basedir,'eli/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/figS2/')
dir.create(savedir,recursive = T)

# load distribution of spearman correlations between trial index and N1/N2/pre stim negative control

data <- readMat(paste0(locations$results_folder,'pub_figs/AllPatientsSpearmanLocsSoz.mat'))
spear.results <- unpack.mat(data,'spear.results')
elec.info <- unpack.mat(data,'elec.info')
ave.cceps <- unpack.mat(data,'ave.cceps')
pts <- name(names(elec.info))

# load SOZ data

soz.all.pts <- unpack.mat.struct(elec.info,'soz')

# rename SOZ to characters
soz.all.pts <- lapply(soz.all.pts, function(x) ifelse(x,yes='SOZ',no='Non-SOZ'))

# tile soz indicator to matrix size so you know if SOZ was stimulated for each CCEP
soz.all.pts.rec <- lapply(soz.all.pts, function(x) matrix(x,nrow=length(x),ncol=length(x))) # check: colSums(soz.all.pts.rec$HUP211=='SOZ') should be 25 for every electrode
soz.all.pts.stim <- lapply(soz.all.pts, function(x) t(matrix(x,nrow=length(x),ncol=length(x)))) # check: rowSums(soz.all.pts.stim$HUP211=='SOZ') should be 25 for every electrode

# load InterIctal spike data

spikes.all.pts <- unpack.mat.struct(elec.info,'spikes')
spikes.all.pts.rec <- lapply(spikes.all.pts, function(x) matrix(x,nrow=length(x),ncol=length(x))) # check: colSums(soz.all.pts.rec$HUP211=='SOZ') should be 25 for every electrode
spikes.all.pts.stim <- lapply(spikes.all.pts, function(x) t(matrix(x,nrow=length(x),ncol=length(x)))) # check: rowSums(soz.all.pts.stim$HUP211=='SOZ') should be 25 for every electrod

# make null spike maps

nperms <- 10
spikes.all.pts.null <- lapply(1:nperms, function(P) 
  lapply(spikes.all.pts, function(spikes.pt) sample(spikes.pt)))
spikes.all.pts.stim.null <- lapply(spikes.all.pts.null, function(spikes.null) 
  lapply(spikes.null, function(x) t(matrix(x,nrow=length(x),ncol=length(x))) ))
spikes.all.pts.rec.null <- lapply(spikes.all.pts.null, function(spikes.null) 
  lapply(spikes.null, function(x) matrix(x,nrow=length(x),ncol=length(x)) ))

# load wm elecs and just use that as a control variable
wm.all.pts <- unpack.mat.struct(elec.info,'wm')

# rename GM to characters
gm.mod.all.pts <- lapply(wm.all.pts, function(x) ifelse(x,yes='WM',no='Non-WM'))

# tile GM indicator to matrix size so you know if GM was stimulated for each CCEP
gm.mod.all.pts.rec <- lapply(gm.mod.all.pts, function(x) matrix(x,nrow=length(x),ncol=length(x))) # check: colSums(gm.mod.all.pts.rec$HUP211=='GM') should be 25 for every electrode
gm.mod.all.pts.stim <- lapply(gm.mod.all.pts, function(x) t(matrix(x,nrow=length(x),ncol=length(x)))) # check: rowSums(gm.mod.all.pts.stim$HUP211=='GM') should be 25 for every electrode

# load phase correlation data
data.hpc.phase <- readMat(paste0(locations$results_folder,'pub_figs/fig4/AllPatientsHippocampalPhaseCLcorr.mat'))
clcorr.results <- unpack.mat(data.hpc.phase,'clcorr.results')

# loop through wave forms, look at number of significant spearman effects for SOZ and Non-SOZ electrodes

waveforms <- c('N1','N2','PreStim')
df.all <- list()
df.null.stim.all <- df.null.rec.all <- list()
wave <- 'N1'

for(wave in waveforms){
  print(wave)
  # unpack mat struct
  data.wave <- unpack.mat(spear.results,wave)
  good.cceps.all.pts <- unpack.mat.struct(data.wave,'good.cceps')
  
  data.wave.clcorr <- unpack.mat(clcorr.results,wave)
  pts.elecs <- lapply(pts, function(pt) names(unpack.mat(data.wave.clcorr,pt)))
  clcorr.all.pts <- unpack.mat.struct.variable(data.wave.clcorr,pts,pts.elecs,'clcorr.trials')
  p.adj.all.pts <- unpack.mat.struct.variable(data.wave.clcorr,pts,pts.elecs,'clcorr.trials.p.adj')
  
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
  
  pt.indicator <- lapply(pts, function(pt) matrix(pt,nrow=nrow(D.all.pts[[pt]]),ncol=ncol(D.all.pts[[pt]])))
  pt.indicator.elecs <- lapply(pts, function(pt) rep(pt,length(chLabels.all[[pt]])))
  
  df <- data.frame(rho=list.mat.to.vec(list.mat.to.vec(clcorr.all.pts)),
                   g=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(good.cceps.all.pts,pts,pts.elecs))),
                   p.adj=list.mat.to.vec(list.mat.to.vec(p.adj.all.pts)),
                   BN.rec=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(Brainnetome.all.pts.rec,pts,pts.elecs))),
                   BN.stim=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(Brainnetome.all.pts.stim,pts,pts.elecs))),
                   D=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(D.all.pts,pts,pts.elecs))),
                   stim.name=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(stim.name,pts,pts.elecs))),
                   stim.name.ana=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(stim.name.ana,pts,pts.elecs))),
                   rec.name=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(rec.name,pts,pts.elecs))),
                   rec.name.ana=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(rec.name.ana,pts,pts.elecs))),
                   pt=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(pt.indicator,pts,pts.elecs))),
                   soz.stim=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(soz.all.pts.stim,pts,pts.elecs))),
                   soz.rec=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(soz.all.pts.rec,pts,pts.elecs))),
                   gm.stim=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(gm.mod.all.pts.stim,pts,pts.elecs))),
                   gm.rec=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(gm.mod.all.pts.rec,pts,pts.elecs))),
                   spikes.stim=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(spikes.all.pts.stim,pts,pts.elecs))),
                   spikes.rec=list.mat.to.vec(list.mat.to.vec(expand.mat.struct(spikes.all.pts.rec,pts,pts.elecs))),
                   wave=wave)
  
  df$soz.non.sz <- ifelse(df$soz.stim == 'Non-SOZ' & df$soz.rec == 'Non-SOZ',yes = 'Outside',no='Inside')
  df$soz.stim.rec <- paste0(df$soz.stim,' -> ',df$soz.rec)
  
  df$Hippocampus.Stim <- df$BN.stim %in% c(215:218) | df$stim.name.ana == 'HIPP'
  df$Hippocampus.Rec <- df$BN.rec %in% c(215:218) | df$rec.name.ana == 'HIPP'
  
  df$Amygdala.Stim <- df$BN.stim %in% c(211:214) | df$stim.name.ana == 'AM'
  df$Amygdala.Rec <- df$BN.rec %in% c(211:214) | df$rec.name.ana == 'AM'
  
  df.good <- df[df$g==1,] # work with smaller df with no NAs to improve speed
  
  df.all[[wave]] <- df
}


# do tests at the CCEP level for each patient to see if CCEPs, rho, failure rate, all localize to SOZ,
# either for quadrants of CCEP matrix or for aggregation into Inside vs. Outside SOZ

df.plt <- do.call(rbind,df.all)
df.plt <- df.plt[df.plt$g == 1 & !is.na(df.plt$rho),]

wave <- 'N1'
pt <- 'HUP223'
results.soz <- list()
df.plt.pt.all <- list()

#waveforms <- c('N1','N2')
coef.names <- c('D','gm.stim','gm.rec','spikes.stim','spikes.rec')
for(wave in waveforms){
  for(pt in pts){
    df.plt.pt <- df.plt[df.plt$pt==pt & df.plt$wave == wave,]
    
    #df.plt.pt <- df.plt.pt[df.plt.pt$rho<0,]
    #df.plt.pt$spikes.rec <- rank_INT(df.plt.pt$spikes.rec)
    #df.plt.pt$spikes.stim <- rank_INT(df.plt.pt$spikes.stim)
    if(nrow(df.plt.pt)>0 & !all(is.na(df.plt.pt$spikes.rec))){
      
      m.rho.sozinout <- lm(rho~log10(D)+gm.stim+gm.rec+soz.non.sz,data=df.plt.pt)
      df.plt.pt$rho.soz.pr <- get.partial.resids.factor(m.rho.sozinout,'soz.non.sz','Outside')$y
      results.soz[[wave]]$rho.soz.inout[[pt]] <- get.coef.p.val(m.rho.sozinout,'soz.non.szOutside')
      
      df.plt.pt.all[[wave]][[pt]] <- df.plt.pt
    }
    
  }
}

lapply(results.soz,as.data.frame)
lapply(results.perm.stim$p,as.data.frame)
lapply(results.perm.rec$p,as.data.frame)

# plot spikes vs. SOZ

df.p.all <- list()
for(test in names(results.soz$N1)){
  df.p <- as.data.frame(sapply(results.soz, function(X) X[[test]]))
  df.p$pt <- rownames(df.p)
  df.p <- collapse.columns(df.p,cnames = waveforms,groupby='pt')
  df.p$values <- p.adjust(df.p$values,method='fdr')
  df.p.all[[test]] <- rename.columns(df.p,c('names','group'),c('wave','pt'))
}

# plot results of CCEP level analysis of SOZ vs. CCEP amplitude (categorical variable, still use partial residuals)
plot.cfg <- list(list(ttl='HippocampalCLCorr',
                      m.name='rho.soz.inout',
                      y.name='rho.soz.pr',
                      x.name='soz.non.sz',
                      xl='SOZ',
                      yl='Rho'))

p.list <- list()
for(cfg in plot.cfg){
  for(wave in waveforms){
    
    df.plt.wave <- do.call(rbind,df.plt.pt.all[[wave]])
    df.p.wave <- df.p.all[[cfg$m.name]][df.p.all[[cfg$m.name]]$wave==wave,]
    df.p.wave$values <- p.signif(df.p.wave$values)
    
    p.list[[wave]][[cfg$ttl]] <- ggplot(df.plt.wave) + theme_classic() + ggtitle(wave)+
      geom_boxplot(aes_string(fill=cfg$x.name,y=cfg$y.name,x='pt'),size=0.1,outlier.size=0.2,outlier.alpha = 0.5,outlier.stroke=0) + standard_plot_addon() +
      geom_text(data=df.p.wave,aes(label=values,y=Inf,x=pt),vjust=1,hjust=1,size=1.5) +
      scale_fill_manual(values=wes_palettes$BottleRocket2[c(2,3)])+
      xlab(cfg$xl) + ylab(cfg$yl)+
      theme(axis.text.x = element_text(angle=90,vjust=0.5)) +standard_plot_addon() +
      theme(plot.margin = margin(0,0,3,0,'mm'))+
      theme(legend.position = 'none',legend.key.size = unit(0.1,'cm'),
            legend.title = element_blank(),
            axis.title.x = element_blank())
    ggsave(plot = p.list[[wave]][[cfg$ttl]],filename = paste0(savedir,cfg$ttl,wave,'.pdf'),width = 4,height=4.5,units= 'cm',useDingbats=FALSE)
    
  }
}

