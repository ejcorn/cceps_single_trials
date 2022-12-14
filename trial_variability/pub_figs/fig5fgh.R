# test for relationship between CCEP amplitude, Rho, Fail rate and amount of spikes in each stimulating or recording
# electrode, using linear model with permutation test to shuffle the electrodes and generate null distribution of betas

rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/miscfxns/statfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig5/fig5fgh/')
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

# loop through wave forms, look at number of significant spearman effects for SOZ and Non-SOZ electrodes

waveforms <- c('N1','N2','PreStim')
df.all <- list()
df.null.stim.all <- df.null.rec.all <- list()
wave <- 'N1'

for(wave in waveforms){
  print(wave)
  # unpack mat struct
  data.wave <- unpack.mat(spear.results,wave)
  spear.all.pts <- unpack.mat.struct(data.wave,'spear.trials')
  good.cceps.all.pts <- unpack.mat.struct(data.wave,'good.cceps')
  p.adj.all.pts <- unpack.mat.struct(data.wave,'spear.trials.p.adj')
  fail.rate <- unpack.mat.struct(data.wave,'fail.rate')
  std.trials <- unpack.mat.struct(data.wave,'std.trials')
  ave.cceps.wave <- unpack.mat.struct(ave.cceps,wave)
  
  if(wave == 'PreStim'){ # make all Nans if prestim because we never calculate an average pre stim network
    ave.cceps.wave <- unpack.mat.struct(ave.cceps,'N1')
    ave.cceps.wave <- lapply(ave.cceps.wave, function(X) matrix(NA,nrow=nrow(X),ncol=ncol(X)))
  }
  
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
  
  pt.indicator <- lapply(pts, function(pt) matrix(pt,nrow=nrow(spear.all.pts[[pt]]),ncol=ncol(spear.all.pts[[pt]])))
  pt.indicator.elecs <- lapply(pts, function(pt) rep(pt,length(chLabels.all[[pt]])))
  
  df <- data.frame(rho=list.mat.to.vec(spear.all.pts),
                   g=list.mat.to.vec(good.cceps.all.pts),
                   p.adj=list.mat.to.vec(p.adj.all.pts),
                   BN.rec=list.mat.to.vec(Brainnetome.all.pts.rec),
                   BN.stim=list.mat.to.vec(Brainnetome.all.pts.stim),
                   D=list.mat.to.vec(D.all.pts),
                   stim.name=list.mat.to.vec(stim.name),
                   stim.name.ana=list.mat.to.vec(stim.name.ana),
                   rec.name=list.mat.to.vec(rec.name),
                   rec.name.ana=list.mat.to.vec(rec.name.ana),
                   pt=list.mat.to.vec(pt.indicator),
                   soz.stim=list.mat.to.vec(soz.all.pts.stim),
                   soz.rec=list.mat.to.vec(soz.all.pts.rec),
                   gm.stim=list.mat.to.vec(gm.mod.all.pts.stim),
                   gm.rec=list.mat.to.vec(gm.mod.all.pts.rec),
                   spikes.stim=list.mat.to.vec(spikes.all.pts.stim),
                   spikes.rec=list.mat.to.vec(spikes.all.pts.rec),
                   ccep=list.mat.to.vec(ave.cceps.wave),
                   fail.rate = list.mat.to.vec(fail.rate),
                   std.trials = list.mat.to.vec(std.trials),
                   wave=wave)
  
  df$soz.non.sz <- ifelse(df$soz.stim == 'Non-SOZ' & df$soz.rec == 'Non-SOZ',yes = 'Outside',no='Inside')
  df$soz.stim.rec <- paste0(df$soz.stim,' -> ',df$soz.rec)
  df$rho <- fisher.r.to.z(df$rho)
  
  df$Hippocampus.Stim <- df$BN.stim %in% c(215:218) | df$stim.name.ana == 'HIPP'
  df$Hippocampus.Rec <- df$BN.rec %in% c(215:218) | df$rec.name.ana == 'HIPP'
  
  df$Amygdala.Stim <- df$BN.stim %in% c(211:214) | df$stim.name.ana == 'AM'
  df$Amygdala.Rec <- df$BN.rec %in% c(211:214) | df$rec.name.ana == 'AM'
  
  # make null
  if(wave != 'PreStim'){
    df.null.stim <- df.null.rec <- df
    for(perm in 1:nperms){
      df.null.stim$spikes.stim <- list.mat.to.vec(spikes.all.pts.stim.null[[perm]])
      df.null.rec$spikes.rec <- list.mat.to.vec(spikes.all.pts.rec.null[[perm]])
      
      df.null.stim.all[[wave]][[perm]] <- df.null.stim[df.null.stim$g==1,]
      df.null.rec.all[[wave]][[perm]] <- df.null.rec[df.null.rec$g==1,]
    }
  }
  
  df.good <- df[df$g==1,] # work with smaller df with no NAs to improve speed
  
  df.all[[wave]] <- df
}


# do tests at the CCEP level for each patient to see if CCEPs, rho, failure rate, all localize to SOZ,
# either for quadrants of CCEP matrix or for aggregation into Inside vs. Outside SOZ

df.plt <- do.call(rbind,df.all)
df.plt <- df.plt[df.plt$g == 1,]

wave <- 'N1'
pt <- 'HUP223'
results.rec <- list()
results.stim <- list()
results.perm.stim <- list()
results.perm.rec <- list()
df.plt.pt.all <- list()

waveforms <- c('N1','N2')
coef.names <- c('D','gm.stim','gm.rec','spikes.stim','spikes.rec')
for(wave in waveforms){
  for(pt in pts){
    df.plt.pt <- df.plt[df.plt$pt==pt & df.plt$wave == wave,]
    df.null.stim.pt <- lapply(df.null.stim.all[[wave]], function(X) X[X$pt==pt,])
    df.null.rec.pt <- lapply(df.null.rec.all[[wave]], function(X) X[X$pt==pt,])
    
    # check
    if(nrow(df.plt.pt) != nrow(df.null.stim.pt[[1]]) | nrow(df.plt.pt) != nrow(df.null.rec.pt[[1]])){stop('null model not aligned with real data')}
    
    #df.plt.pt <- df.plt.pt[df.plt.pt$rho<0,]
    #df.plt.pt$spikes.rec <- rank_INT(df.plt.pt$spikes.rec)
    #df.plt.pt$spikes.stim <- rank_INT(df.plt.pt$spikes.stim)
    if(nrow(df.plt.pt)>0 & !all(is.na(df.plt.pt$spikes.rec))){
      
      m.rho.spikes.stim <- lm(spikes.stim~log10(D)+gm.stim+gm.rec+rho ,data=df.plt.pt)

      results.stim$p[[wave]]$rho.spikes[[pt]] <- get.coef.p.val(m.rho.spikes.stim,'rho')
      results.stim$b[[wave]]$rho.spikes[[pt]] <- get.coef.beta(m.rho.spikes.stim,'rho')
      results.stim$m[[wave]]$rho.spikes[[pt]] <- m.rho.spikes.stim
      results.perm.stim$p[[wave]]$rho.spikes[[pt]] <- lm.perm.null.df.list(m.rho.spikes.stim,df.null.stim.pt,coef.names)['rho']
      
      m.rho.spikes.rec <- lm(spikes.rec~log10(D)+gm.stim+gm.rec+rho,data=df.plt.pt)
      results.rec$p[[wave]]$rho.spikes[[pt]] <- get.coef.p.val(m.rho.spikes.rec,'rho')
      results.rec$b[[wave]]$rho.spikes[[pt]] <- get.coef.beta(m.rho.spikes.rec,'rho')
      results.rec$m[[wave]]$rho.spikes[[pt]] <- m.rho.spikes.rec
      results.perm.rec$p[[wave]]$rho.spikes[[pt]] <- lm.perm.null.df.list(m.rho.spikes.rec,df.null.rec.pt,coef.names)['rho']
      
      m.ccep.spikes.rec <- lm(spikes.rec~log10(D)+gm.stim+gm.rec+log10(ccep),data=df.plt.pt)
      results.rec$p[[wave]]$ccep.spikes[[pt]] <- get.coef.p.val(m.ccep.spikes.rec,'log10(ccep)')
      results.rec$b[[wave]]$ccep.spikes[[pt]] <- get.coef.beta(m.ccep.spikes.rec,'log10(ccep)')
      results.rec$m[[wave]]$ccep.spikes[[pt]] <- m.ccep.spikes.rec
      results.perm.rec$p[[wave]]$ccep.spikes[[pt]] <- lm.perm.null.df.list(m.ccep.spikes.rec,df.null.rec.pt,coef.names)['log10(ccep)']
      
      m.ccep.spikes.stim <- lm(spikes.stim~log10(D)+gm.stim+gm.rec+log10(ccep),data=df.plt.pt)
      results.stim$p[[wave]]$ccep.spikes[[pt]] <- get.coef.p.val(m.ccep.spikes.stim,'log10(ccep)')
      results.stim$b[[wave]]$ccep.spikes[[pt]] <- get.coef.beta(m.ccep.spikes.stim,'log10(ccep)')
      results.stim$m[[wave]]$ccep.spikes[[pt]] <- m.ccep.spikes.stim
      results.perm.stim$p[[wave]]$ccep.spikes[[pt]] <- lm.perm.null.df.list(m.ccep.spikes.stim,df.null.stim.pt,coef.names)['log10(ccep)']
      
      m.std.spikes.stim <- lm(spikes.stim~log10(std.trials)+log10(D)+gm.stim+gm.rec,data=df.plt.pt)
      results.stim$p[[wave]]$std.spikes[[pt]] <- get.coef.p.val(m.std.spikes.stim,'log10(std.trials)')
      results.stim$b[[wave]]$std.spikes[[pt]] <- get.coef.beta(m.std.spikes.stim,'log10(std.trials)')
      results.stim$m[[wave]]$std.spikes[[pt]] <- m.std.spikes.stim
      results.perm.stim$p[[wave]]$std.spikes[[pt]] <- lm.perm.null.df.list(m.std.spikes.stim,df.null.stim.pt,coef.names)['log10(std.trials)']

      m.std.spikes.rec<- lm(spikes.rec~log10(std.trials)+log10(D)+gm.stim+gm.rec,data=df.plt.pt)      
      results.rec$p[[wave]]$std.spikes[[pt]] <- get.coef.p.val(m.std.spikes.rec,'log10(std.trials)')
      results.rec$b[[wave]]$std.spikes[[pt]] <- get.coef.beta(m.std.spikes.rec,'log10(std.trials)')
      results.rec$m[[wave]]$std.spikes[[pt]] <- m.std.spikes.rec
      results.perm.rec$p[[wave]]$std.spikes[[pt]] <- lm.perm.null.df.list(m.std.spikes.rec,df.null.rec.pt,coef.names)['log10(std.trials)']
      
      # +Hippocampus.Stim+Amygdala.Stim+Hippocampus.Rec+Amygdala.Rec
      df.plt.pt.all[[wave]][[pt]] <- df.plt.pt
    } else{
   
      
    }
    
  }
  
}

lapply(results.stim$p,as.data.frame)
lapply(results.perm.stim$p,as.data.frame)
lapply(results.perm.rec$p,as.data.frame)

lapply(results.stim$m$N1, function(X) 
  sapply(X, function(mdls) summary(mdls)$r.sq))

# plot spikes vs. SOZ

df.p.stim.all <- df.p.rec.all <- list()
df.b.stim.all <- df.b.rec.all <- list()
for(test in names(results.stim$p$N1)){
  df.p <- as.data.frame(sapply(results.stim$p, function(X) X[[test]]))
  df.p$pt <- rownames(df.p)
  df.p <- collapse.columns(df.p,cnames = waveforms,groupby='pt')
  df.p$stim.rec <- 'stim'
  df.p$values <- p.adjust(df.p$values,method='fdr')
  df.p.stim.all[[test]] <- rename.columns(df.p,c('names','group'),c('wave','pt'))
  
  df.p <- as.data.frame(sapply(results.rec$p, function(X) X[[test]]))
  df.p$pt <- rownames(df.p)
  df.p <- collapse.columns(df.p,cnames = waveforms,groupby='pt')
  df.p$stim.rec <- 'rec'
  df.p$values <- p.adjust(df.p$values,method='fdr')
  df.p.rec.all[[test]] <- rename.columns(df.p,c('names','group'),c('wave','pt'))
  
  df.b <- as.data.frame(sapply(results.stim$b, function(X) X[[test]]))
  df.b$pt <- rownames(df.b)
  df.b <- collapse.columns(df.b,cnames = waveforms,groupby='pt')
  df.b$stim.rec <- 'stim'
  df.b.stim.all[[test]] <- rename.columns(df.b,c('names','group'),c('wave','pt'))
  
  df.b <- as.data.frame(sapply(results.rec$b, function(X) X[[test]]))
  df.b$pt <- rownames(df.b)
  df.b <- collapse.columns(df.b,cnames = waveforms,groupby='pt')
  df.b$stim.rec <- 'rec'
  df.b.rec.all[[test]] <- rename.columns(df.b,c('names','group'),c('wave','pt'))
  
}

# plot results of CCEP level analysis of spike rate vs. rho
plot.cfg <- list(list(ttl='RhoVsSpikes',
                      p.ttl= 'Rho',
                      m.name='rho.spikes',
                      y.name='rho.spikes.pr',
                      x.name='spikes.rec'),
                 list(ttl='CCEPVsSpikes',
                      p.ttl= 'Amp.',
                      m.name='ccep.spikes',
                      y.name='ccep.spikes.pr',
                      x.name='spikes.rec'),
                 list(ttl='StdVsSpikes',
                      m.name='std.spikes',
                      p.ttl= 'Std. Dev.')
)

p <- list()
for(cfg in plot.cfg){
  for(wave in waveforms){
    
    df.plt.wave <- do.call(rbind,df.plt.pt.all[[wave]])
    df.p.stim.wave <- df.p.stim.all[[cfg$m.name]][df.p.stim.all[[cfg$m.name]]$wave==wave,]
    df.b.stim.wave <- df.b.stim.all[[cfg$m.name]][df.b.stim.all[[cfg$m.name]]$wave==wave,]
    df.p.rec.wave <- df.p.rec.all[[cfg$m.name]][df.p.rec.all[[cfg$m.name]]$wave==wave,]
    df.b.rec.wave <- df.b.rec.all[[cfg$m.name]][df.b.rec.all[[cfg$m.name]]$wave==wave,]
    
    df.p.rec.wave$values <- p.signif(df.p.rec.wave$values,ns='')
    df.p.stim.wave$values <- p.signif(df.p.stim.wave$values,ns='')
    
    df.b.wave <- rbind(df.b.rec.wave,df.b.stim.wave)
    ylim <- 1.15*max(abs(df.b.wave$values),na.rm = T)
    
    p[[wave]][[cfg$ttl]] <- ggplot(df.b.wave,aes(x=pt,y=values,fill=stim.rec)) +  theme_classic() + 
      geom_col(position=position_dodge()) + 
      geom_text(data=df.p.rec.wave,aes(label=values,y=Inf,x=pt),color=wes_palettes$Darjeeling1[2],vjust=1,size=2) +
      geom_text(data=df.p.stim.wave,aes(label=values,y=-Inf,x=pt),color=wes_palettes$Darjeeling1[4],vjust=-1,size=2) +
      scale_fill_manual(values=wes_palettes$Darjeeling1[c(2,4)])+
      scale_y_continuous(limits=c(-ylim,ylim)) +
      ylab(expression(beta))+
      ggtitle(paste0(wave,': ',cfg$p.ttl)) +
      standard_plot_addon() +
      theme(plot.margin= margin(0,0,3,0, "mm"))+
      theme(axis.text.x = element_text(angle=90,vjust=0.5)) +
      theme(legend.position = 'none',
            legend.key.size = unit(0.1,'cm'),
            legend.title = element_blank(),
            axis.title.x = element_blank())
    ggsave(plot = p[[wave]][[cfg$ttl]],filename = paste0(savedir,cfg$ttl,wave,'.pdf'),width = 3,height=4,units= 'cm',useDingbats=FALSE)
    
  }
}

p.all <- plot_grid(plotlist=list(p$N1$CCEPVsSpikes,p$N2$CCEPVsSpikes,
                                 p$N1$RhoVsSpikes,p$N2$RhoVsSpikes,
                                 p$N1$StdVsSpikes,p$N2$StdVsSpikes),
                   align='hv',ncol=2)
ggsave(plot = p.all,filename = paste0(savedir,'Fig5fgh.pdf'),width = 8.75,height=10,units= 'cm',useDingbats=FALSE)
