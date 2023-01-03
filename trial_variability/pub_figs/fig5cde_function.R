# this script makes a glm for each patient to predict whether each CCEP belongs to the SOZ based 
# on inter-CCEP distance, GM/WM, CCEP amplitude and standard deviation of CCEP amplitude
# combines N1 and N2 into the same model
# then you will need to add 

rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/miscfxns/statfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig5/fig5cde/')
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

# make null SOZ maps

nperms <- 10
soz.all.pts.null <- lapply(1:nperms, function(P) 
  lapply(soz.all.pts, function(soz.pt) sample(soz.pt)))
soz.all.pts.stim.null <- lapply(soz.all.pts.null, function(soz.null) 
  lapply(soz.null, function(x) t(matrix(x,nrow=length(x),ncol=length(x))) ))
soz.all.pts.rec.null <- lapply(soz.all.pts.null, function(soz.null) 
  lapply(soz.null, function(x) matrix(x,nrow=length(x),ncol=length(x)) ))

soz.non.sz.all.pts.null <- lapply(1:nperms, function(x) NULL)
for(P in 1:nperms){
  for(pt in pts){
    soz.non.sz.all.pts.null[[P]][[pt]] <- ifelse(soz.all.pts.stim.null[[P]][[pt]] == 'Non-SOZ' & soz.all.pts.rec.null[[P]][[pt]] == 'Non-SOZ',yes = 'Outside',no='Inside')
  }
}

rm(list=c('soz.all.pts.stim.null','soz.all.pts.rec.null'))
# load InterIctal spike data

spikes.all.pts <- unpack.mat.struct(elec.info,'spikes')
spikes.all.pts.rec <- lapply(spikes.all.pts, function(x) matrix(x,nrow=length(x),ncol=length(x))) # check: colSums(soz.all.pts.rec$HUP211=='SOZ') should be 25 for every electrode
spikes.all.pts.stim <- lapply(spikes.all.pts, function(x) t(matrix(x,nrow=length(x),ncol=length(x)))) # check: rowSums(soz.all.pts.stim$HUP211=='SOZ') should be 25 for every electrod

# load wm elecs and just use that as a control variable
wm.all.pts <- unpack.mat.struct(elec.info,'wm')

# rename GM to characters
gm.mod.all.pts <- lapply(wm.all.pts, function(x) ifelse(x,yes='WM',no='Non-WM'))

# tile GM indicator to matrix size so you know if GM was stimulated for each CCEP
gm.mod.all.pts.rec <- lapply(gm.mod.all.pts, function(x) matrix(x,nrow=length(x),ncol=length(x))) # check: colSums(gm.mod.all.pts.rec$HUP211=='GM') should be 25 for every electrode
gm.mod.all.pts.stim <- lapply(gm.mod.all.pts, function(x) t(matrix(x,nrow=length(x),ncol=length(x)))) # check: rowSums(gm.mod.all.pts.stim$HUP211=='GM') should be 25 for every electrodee

# loop through wave forms, look at number of significant spearman effects for SOZ and Non-SOZ electrodes

waveforms <- c('N1','N2','PreStim')
df.all <- list()
df.null.all <- list()
results <- list()
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
  
  df <- data.frame(rho=fisher.r.to.z(list.mat.to.vec(spear.all.pts)),
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
  
  df$Hippocampus.Stim <- df$BN.stim %in% c(215:218) | df$stim.name.ana == 'Hippocampus'
  df$Hippocampus.Rec <- df$BN.rec %in% c(215:218) | df$rec.name.ana == 'Hippocampus'
  
  df$Amygdala.Stim <- df$BN.stim %in% c(211:214) | df$stim.name.ana == 'Amygdala'
  df$Amygdala.Rec <- df$BN.rec %in% c(211:214) | df$rec.name.ana == 'Amygdala'
  
  df$Motor.Stim <- df$BN.stim %in% c(53:64) | df$stim.name.ana == 'precentral'
  df$Motor.Rec <- df$BN.rec %in% c(53:64) | df$rec.name.ana == 'precentral'
  
  df$Somatosensory.Stim <- df$BN.stim %in% c(155:162) | df$stim.name.ana == 'postcentral'
  df$Somatosensory.Rec <- df$BN.rec %in% c(155:162) | df$rec.name.ana == 'postcentral'  
  
  df$Visual.Stim <- df$BN.stim %in% c(189:210) | df$stim.name.ana == 'occipital'
  df$Visual.Rec <- df$BN.rec %in% c(189:210) | df$rec.name.ana == 'occipital'    
  
  # make null
  df.null <- df
  for(perm in 1:nperms){
    df.null$soz.non.sz <- list.mat.to.vec(soz.non.sz.all.pts.null[[perm]])
    df.null.all[[wave]][[perm]] <- df.null[df.null$g==1 & !is.na(df.null$soz.stim),]
  }
  
  df.good <- df[df$g==1,] # work with smaller df with no NAs to improve speed
  
  df.all[[wave]] <- df
}


# do tests at the CCEP level for each patient to see if CCEPs, rho, failure rate, all localize to SOZ,
# either for quadrants of CCEP matrix or for aggregation into Inside vs. Outside SOZ

df.plt <- do.call(rbind,df.all)
df.plt <- df.plt[df.plt$g == 1,]

sapply(df.all,function(X) sum(X$g==1))
wave <- 'N1'
waveforms <- c('N1','N2')
pt <- 'HUP223'
results.soz <- list()
results.soz.perm <- list()
m.soz.all <- list()
m.soz.pred <- list()
m.soz.pred.auc <- list()
df.plt.pt.all <- list()
coef.names <- c('D','gm.stim','gm.rec','soz.non.sz')

for(wave in waveforms){
  for(pt in pts){
    df.plt.pt <- df.plt[df.plt$pt==pt & !is.na(df.plt$soz.stim) & df.plt$wave == wave,]
    df.null.pt <- lapply(df.null.all[[wave]], function(X) X[X$pt==pt,])
    
    if(nrow(df.plt.pt) != nrow(df.null.pt[[1]])){stop('null model not aligned with real data')}
    
    #df.plt.pt <- df.plt.pt[df.plt.pt$rho<0,]
    #df.plt.pt$spikes.rec <- rank_INT(df.plt.pt$spikes.rec)
    #df.plt.pt$spikes.stim <- rank_INT(df.plt.pt$spikes.stim)
    if(nrow(df.plt.pt)>0){
      if(length(unique(df.plt.pt$soz.non.sz[!is.na(df.plt.pt$ccep)]))>1){
        
        m.ccep.sozinout <- lm(log10(ccep)~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+soz.non.sz,data=df.plt.pt)
        df.plt.pt$ccep.soz.pr <- get.partial.resids.factor(m.ccep.sozinout,'soz.non.sz','Outside')$y
        m.soz.all[[wave]]$ccep[[pt]] <- m.ccep.sozinout
        
        if(wave != 'PreStim'){
          
          results.soz[[wave]]$ccep.soz.inout[[pt]] <- get.coef.p.val(m.ccep.sozinout,'soz.non.szOutside')
          results.soz.perm[[wave]]$ccep.soz.inout[[pt]] <- lm.perm.null.df.list(m.ccep.sozinout,df.null.pt,coef.names)['soz.non.szOutside']
          
        } else{
          results.soz[[wave]]$ccep.soz.inout[[pt]] <- results.soz.perm[[wave]]$ccep.soz.inout[[pt]] <- NA
        }
        
      } else{
        df.plt.pt$ccep.soz.pr <- NA
        results.soz[[wave]]$ccep.soz.inout[[pt]] <- results.soz.perm[[wave]]$ccep.soz.inout[[pt]] <- NA
        results.soz[[wave]]$rho.soz.inout[[pt]] <- results.soz.perm[[wave]]$rho.soz.inout[[pt]] <- NA
        results.soz[[wave]]$fail.rate.soz.inout[[pt]] <- results.soz.perm[[wave]]$fail.rate.soz.inout[[pt]] <- NA
      }
      
      m.rho.sozinout <- lm(rho~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+soz.non.sz,data=df.plt.pt)
      df.plt.pt$rho.soz.pr <- get.partial.resids.factor(m.rho.sozinout,'soz.non.sz','Outside')$y
      m.soz.all[[wave]]$rho[[pt]] <- m.rho.sozinout
      results.soz[[wave]]$rho.soz.inout[[pt]] <- get.coef.p.val(m.rho.sozinout,'soz.non.szOutside')
      results.soz.perm[[wave]]$rho.soz.inout[[pt]] <- lm.perm.null.df.list(m.rho.sozinout,df.null.pt,coef.names)['soz.non.szOutside']

      
      m.std.trials.sozinout <- lm(std.trials~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+soz.non.sz,data=df.plt.pt)
      df.plt.pt$std.trials.soz.pr <- get.partial.resids.factor(m.std.trials.sozinout,'soz.non.sz','Outside')$y
      m.soz.all[[wave]]$std[[pt]] <- m.std.trials.sozinout
      results.soz[[wave]]$std.trials.soz.inout[[pt]] <- get.coef.p.val(m.std.trials.sozinout,'soz.non.szOutside')
      results.soz.perm[[wave]]$std.trials.soz.inout[[pt]] <- lm.perm.null.df.list(m.std.trials.sozinout,df.null.pt,coef.names)['soz.non.szOutside']
      
      m.log.std.trials.sozinout <- lm(log10(std.trials)~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+soz.non.sz,data=df.plt.pt)
      df.plt.pt$log.std.trials.soz.pr <- get.partial.resids.factor(m.log.std.trials.sozinout,'soz.non.sz','Outside')$y
      m.soz.all$log.std[[pt]] <- m.log.std.trials.sozinout
      results.soz[[wave]]$log.std.trials.soz.inout[[pt]] <- get.coef.p.val(m.log.std.trials.sozinout,'soz.non.szOutside')
      results.soz.perm[[wave]]$log.std.trials.soz.inout[[pt]] <- lm.perm.null.df.list(m.log.std.trials.sozinout,df.null.pt,coef.names)['soz.non.szOutside']
      
      df.plt.pt$soz.non.sz.binary <- df.plt.pt$soz.non.sz == 'Inside'
      m.soz.pred[[wave]]$std[[pt]] <- glm(soz.non.sz.binary~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+scale(std.trials),data=df.plt.pt)
      m.soz.pred[[wave]]$rho[[pt]] <- glm(soz.non.sz.binary~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+scale(rho),data=df.plt.pt)
      m.soz.pred[[wave]]$ccep[[pt]] <- glm(soz.non.sz.binary~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+scale(log10(ccep)),data=df.plt.pt)
      #m.soz.pred.auc[[wave]][[pt]] <- roc(response=df.plt.pt$soz.non.sz.binary,predictor = predict(m.soz.pred[[wave]][[pt]],df.plt.pt,type='response'))
    }
    df.plt.pt.all[[wave]][[pt]] <- df.plt.pt
  }
  
}

save(m.soz.pred,file=paste0(savedir,'Fig5cdeSOZOutcomeModels.RData'))

lapply(results.soz,as.data.frame)
lapply(results.soz.perm,as.data.frame)

df.p.all <- list()
for(test in names(results.soz$N1)){
  df.p <- as.data.frame(sapply(results.soz, function(X) X[[test]]))
  df.p$pt <- rownames(df.p)
  df.p <- collapse.columns(df.p,cnames = waveforms,groupby='pt')
  df.p$values <- p.adjust(df.p$values,method='fdr')
  df.p.all[[test]] <- rename.columns(df.p,c('names','group'),c('wave','pt'))
}

# plot results of CCEP level analysis of SOZ vs. CCEP amplitude (categorical variable, still use partial residuals)
plot.cfg <- list(list(ttl='CCEPVsSOZ',
                      m.name='ccep.soz.inout',
                      y.name='ccep.soz.pr',
                      x.name='soz.non.sz',
                      xl='SOZ',
                      yl='log10(Amp.)'),
                 list(ttl='RhoVsSOZ',
                      m.name='rho.soz.inout',
                      y.name='rho.soz.pr',
                      x.name='soz.non.sz',
                      xl='SOZ',
                      yl='Rho'),
                 list(ttl='StdVsSOZ',
                      m.name='std.trials.soz.inout',
                      y.name='std.trials.soz.pr',
                      x.name='soz.non.sz',
                      xl='SOZ',
                      yl='Std. Dev.'))

p.list <- list()
for(cfg in plot.cfg){
  for(wave in waveforms){
    
    df.plt.wave <- do.call(rbind,df.plt.pt.all[[wave]])
    df.p.wave <- df.p.all[[cfg$m.name]][df.p.all[[cfg$m.name]]$wave==wave,]
    df.p.wave$values <- p.signif(df.p.wave$values,ns='')
    
    p.list[[wave]][[cfg$ttl]] <- ggplot(df.plt.wave) + theme_classic() + ggtitle(wave)+
      geom_boxplot(aes_string(fill=cfg$x.name,y=cfg$y.name,x='pt'),size=0.1,outlier.size=0.2,outlier.alpha = 0.5,outlier.stroke=0) + standard_plot_addon() +
      geom_text(data=df.p.wave,aes(label=values,y=Inf,x=pt),vjust=1,hjust=1,size=2) +
      scale_fill_manual(values=wes_palettes$BottleRocket2[c(2,3)])+
      xlab(cfg$xl) + ylab(cfg$yl)+
      theme(axis.text.x = element_text(angle=90,vjust=0.5)) +standard_plot_addon() +
      theme(plot.margin = margin(0,0,3,0,'mm'))+
      theme(legend.position = 'none',legend.key.size = unit(0.1,'cm'),
            legend.title = element_blank(),
            axis.title.x = element_blank())
    ggsave(plot = p.list[[wave]][[cfg$ttl]],filename = paste0(savedir,cfg$ttl,wave,'.pdf'),width = 3,height=4,units= 'cm',useDingbats=FALSE)
    
  }
}

p.all <- plot_grid(plotlist=list(p.list$N1$CCEPVsSOZ,p.list$N2$CCEPVsSOZ,
                                 p.list$N1$RhoVsSOZ,p.list$N2$RhoVsSOZ,
                                 p.list$N1$StdVsSOZ,p.list$N2$StdVsSOZ),
                   align='hv',ncol=2)
ggsave(plot = p.all,filename = paste0(savedir,'Fig5cdeFunctionControl.pdf'),width = 8.75,height=10,units= 'cm',useDingbats=FALSE)


# Now also show the relationship between function and CCEP metrics

plot.cfg <- list(list(ttl='CCEP',
                      m.name='ccep',
                      y.name='ccep.soz.pr',
                      x.name='soz.non.sz',
                      xl='SOZ',
                      yl='log10(Amp.)'),
                 list(ttl='Rho',
                      m.name='rho',
                      y.name='rho.soz.pr',
                      x.name='soz.non.sz',
                      xl='SOZ',
                      yl='Rho'),
                 list(ttl='Standard Deviation',
                      m.name='std',
                      y.name='std.trials.soz.pr',
                      x.name='soz.non.sz',
                      xl='SOZ',
                      yl='Std. Dev.'))

nice.coef.names <- c(`log10(D)`='log10(D)',gm.stimWM="GM Stim",gm.recWM='GM Rec',Somatosensory.StimTRUE='SOM Stim',
                     Somatosensory.RecTRUE='SOM Rec',Motor.StimTRUE='MOT Stim',Motor.RecTRUE='MOT Rec',
                     Visual.StimTRUE='VIS Stim',Visual.RecTRUE='VIS Rec',
                     soz.non.szOutside='Outside SOZ')
p.list <- list()
for(wave in waveforms){
  for(cfg in plot.cfg){
    # calculate standardized betas for plot
    std.beta <- lapply(m.soz.all[[wave]][[cfg$m.name]],lm.beta)
    # extract betas and p-values for all subjects and all coefficients in one matrix
    b.matrix <- sapply(std.beta,get.lm.std.betas)[-1,]
    p.matrix <- sapply(std.beta,get.lm.p.vals)[-1,]
    # give coefficients nice names for plot
    rownames(b.matrix) <- unname(nice.coef.names[rownames(b.matrix)])
    rownames(p.matrix) <- unname(nice.coef.names[rownames(p.matrix)])
    # correct p values for each coefficient across patient
    p.adj.matrix <- t(sapply(rownames(p.matrix),function(var) p.adjust(p.matrix[var,],method='fdr')))
    # make p-values into * code for significance level
    p.matrix <- p.signif.matrix(p.adj.matrix,ns=' ')
    # make heatmap plot of betas with overlayed asterisk
    p.list[[wave]][[cfg$m.name]] <- imagesc(X=b.matrix,ttl = paste0(wave,': ',cfg$ttl),
            overlay = p.matrix, # rotate text
            overlay.text.col = 'white',
            overlay.text.sz = 2.5,
            overlay.text.angle = 90,
            overlay.text.vjust = 1,
            cmap = 'redblue_asymmetric_nagrey',
            caxis_name = expression(beta),
            clim=c(-0.5,0.2)) + 
      standard_plot_addon() + x_text_90() + nice_cbar()
  }
  p.all <- plot_grid(plotlist=list(p.list[[wave]]$ccep,
                                   p.list[[wave]]$rho,
                                   p.list[[wave]]$std),
                     align='hv',ncol=3)
  ggsave(plot = p.all,filename = paste0(savedir,'Fig5cdeBetas',wave,'.pdf'),width = 18,height=6,units= 'cm',useDingbats=FALSE)
}

##############################
### Make regression tables ###
##############################

pts <- names(m.soz.all[[wave]][[cfg$m.name]])
for(wave in waveforms){
  reg.tab.rho <- reg.tab.ccep <- reg.tab.std <- list()
  p.matrix <- sapply(m.soz.all[[wave]][['rho']],get.lm.p.vals)
  p.adj.matrix.rho <- t(sapply(rownames(p.matrix),function(var) p.adjust(p.matrix[var,],method='fdr')))
  p.matrix <- sapply(m.soz.all[[wave]][['ccep']],get.lm.p.vals)
  p.adj.matrix.ccep <- t(sapply(rownames(p.matrix),function(var) p.adjust(p.matrix[var,],method='fdr')))
  p.matrix <- sapply(m.soz.all[[wave]][['std']],get.lm.p.vals)
  p.adj.matrix.std <- t(sapply(rownames(p.matrix),function(var) p.adjust(p.matrix[var,],method='fdr')))
  
  for(pt in pts){

      m.ccep.sozinout <- lm(log10(ccep)~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+soz.non.sz,data=df.plt.pt)
      m.rho.sozinout <- lm(rho~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+soz.non.sz,data=df.plt.pt)
      m.std.trials.sozinout <- lm(std.trials~log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+soz.non.sz,data=df.plt.pt)

      reg.tab.rho[[pt]] <- assemble.reg.tab(m.soz.all[[wave]][['rho']][[pt]],'Rho',pt,wave,p.adj.matrix.rho)
      reg.tab.rho[[pt]] <- cbind(reg.tab.rho[[pt]],'') # add blank row to horizontally concatenate
      reg.tab.ccep[[pt]] <- assemble.reg.tab(m.soz.all[[wave]][['ccep']][[pt]],'CCEP Amplitude',pt,wave,p.adj.matrix.ccep)
      reg.tab.ccep[[pt]] <- cbind(reg.tab.ccep[[pt]],'') # add blank row to horizontally concatenate
      reg.tab.std[[pt]] <- assemble.reg.tab(m.soz.all[[wave]][['std']][[pt]],'CCEP Amplitude Standard Deviation',pt,wave,p.adj.matrix.std)
      reg.tab.std[[pt]] <- cbind(reg.tab.std[[pt]],'') # add blank row to horizontally concatenate
      
    }
  
  reg.tab.rho <- do.call(cbind,reg.tab.rho)
  write.table(reg.tab.rho,file = paste0(savedir,'RegressionTable_',wave,'RhoSOZ_Fig5abc.csv'),row.names = F,col.names = F,sep=',')
  reg.tab.ccep <- do.call(cbind,reg.tab.ccep)
  write.table(reg.tab.ccep,file = paste0(savedir,'RegressionTable_',wave,'CCEPAmplitudeSOZ_Fig5abc.csv'),row.names = F,col.names = F,sep=',')
  reg.tab.std <- do.call(cbind,reg.tab.std)
  write.table(reg.tab.std,file = paste0(savedir,'RegressionTable_',wave,'CCEPStandardDeviationSOZ_Fig5abc.csv'),row.names = F,col.names = F,sep=',')
  
}
