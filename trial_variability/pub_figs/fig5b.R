# this script is same as Fig5d but controls for anatomic location, GM/WM status of electrodes and inter electrode distance 
# in assessing relationship between CCEP amplitude and SOZ

rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/miscfxns/statfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig5/fig5b/')
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
wave <- 'N1'

for(wave in waveforms){
  print(wave)
  # unpack mat struct
  data.wave <- unpack.mat(spear.results,wave)
  spear.all.pts <- unpack.mat.struct(data.wave,'spear.trials')
  good.cceps.all.pts <- unpack.mat.struct(data.wave,'good.cceps')
  p.adj.all.pts <- unpack.mat.struct(data.wave,'spear.trials.p.adj')
  fail.rate <- unpack.mat.struct(data.wave,'fail.rate')
  coefvar.trials <- unpack.mat.struct(data.wave,'coefvar.trials')
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
                   coefvar = list.mat.to.vec(coefvar.trials),
                   std = list.mat.to.vec(std.trials),
                   wave=wave,
                   stringsAsFactors = F)
  
  df$soz.non.sz <- ifelse(df$soz.stim == 'Non-SOZ' & df$soz.rec == 'Non-SOZ',yes = 'Outside',no='Inside')
  df$soz.stim.rec <- paste0(df$soz.stim,' -> ',df$soz.rec)
  
  df$Hippocampus.Stim <- df$BN.stim %in% c(215:218) | df$stim.name.ana == 'HIPP'
  df$Hippocampus.Rec <- df$BN.rec %in% c(215:218) | df$rec.name.ana == 'HIPP'
  
  df$Amygdala.Stim <- df$BN.stim %in% c(211:214) | df$stim.name.ana == 'AM'
  df$Amygdala.Rec <- df$BN.rec %in% c(211:214) | df$rec.name.ana == 'AM'
  
  df.good <- df[df$g==1,] # work with smaller df with no NAs to improve speed
  
  df.all[[wave]] <- df
}


# count non-artifactual CCEPs we measured - you can't analyze a CCEP where either N1 or N2 is 0 (b/c findpeaks couldn't find peak)
# but, the total number of nonartifactual CCEPs is greater than the CCEPs we analyzed (almost double) 
count.nonartifact.cceps <- sum(df.all$N1$ccep >0 | df.all$N2$ccep>0,na.rm=T)
count.analyzed.cceps <- sum(df.all$N1$ccep >0 & df.all$N2$ccep>0,na.rm=T)

Intra.SOZ.CCEPs <- lapply(df.all, function(X) X[X$soz.stim.rec == 'SOZ -> SOZ' & !is.na(X$ccep),])
Intra.SOZ.CCEPs <- lapply(df.all, function(X) X[X$soz.stim.rec == 'SOZ -> SOZ' & !is.na(X$ccep) & X$ccep!=0,])
lapply(df.all, function(X) sum(X$g[X$soz.stim.rec == 'SOZ -> SOZ']))

# do tests at the CCEP level for each patient to see if CCEPs, rho, failure rate, all localize to SOZ,
# either for quadrants of CCEP matrix or for aggregation into Inside vs. Outside SOZ

df.plt <- do.call(rbind,df.all)
df.plt <- df.plt[df.plt$g == 1,]

wave <- 'N1'
pt <- 'HUP223'
results <- list()
results.lme <- list()
df.plt.pt.all <- list()

coef.names <- c('std','D','ccep','gm.stim','gm.rec')
waveforms <- c('N1','N2')
for(wave in waveforms){

  # results.lme[[wave]]$lme4.std <- lmer(log10(ccep)~log10(D)+std+gm.stim+gm.rec + 1|pt,data=df.plt[df.plt$wave==wave,],na.action = na.omit,REML=T)
  results.lme[[wave]]$lme.std <- lme(log10(ccep)~log10(D)+std+gm.stim+gm.rec,random=~1|pt,data=df.plt[df.plt$wave==wave,],na.action = na.omit,method = 'ML')
  # results.lme[[wave]]$lme4.std <- lmer(log10(ccep)~log10(D)+std+gm.stim+gm.rec + 1|pt,data=df.plt[df.plt$wave==wave,],na.action = na.omit,REML=T)
  # results.lme[[wave]]$lme.rho <- lme(log10(ccep)~log10(D)+rho+gm.stim+gm.rec,random=~1|pt,data=df.plt[df.plt$wave==wave,],na.action = na.omit)
  # 
  data("hirose")
  for(pt in pts){
    df.plt.pt <- df.plt[df.plt$pt==pt & df.plt$wave == wave & !is.na(df.plt$D),]
    df.plt.pt$rho <- fisher.r.to.z(df.plt.pt$rho)

    if(nrow(df.plt.pt)>0){
      
        m.ccep.coefvar <- lm(log10(ccep)~coefvar+log10(D)+gm.stim+gm.rec,data=df.plt.pt)
        results$p[[wave]]$ccep.coefvar[[pt]] <- get.coef.p.val(m.ccep.coefvar,'coefvar')
        results$b[[wave]]$ccep.coefvar[[pt]] <- get.coef.beta(m.ccep.coefvar,'coefvar')
        df.plt.pt$ccep.coefvar.pr <- get.partial.resids(m.ccep.coefvar,'coefvar')$y
        
        m.ccep.std <- lm(log10(ccep)~std+log10(D)+gm.stim+gm.rec,data=df.plt.pt)
        results$p[[wave]]$ccep.std[[pt]] <- get.coef.p.val(m.ccep.std,'std')
        results$b[[wave]]$ccep.std[[pt]] <- get.coef.beta(m.ccep.std,'std')
        df.plt.pt$ccep.std.pr <- get.partial.resids(m.ccep.std,'std')$y

        m.ccep.rho <- lm(log10(ccep)~rho+log10(D)+gm.stim+gm.rec,data=df.plt.pt)
        results$p[[wave]]$ccep.rho[[pt]] <- get.coef.p.val(m.ccep.rho,'rho')
        results$b[[wave]]$ccep.rho[[pt]] <- get.coef.beta(m.ccep.rho,'rho')
        df.plt.pt$ccep.rho.pr <- get.partial.resids(m.ccep.rho,'rho')$y

        if(mean(is.na(df.plt.pt$fail.rate) | df.plt.pt$fail.rate == 0) <0.5){
          
          m.ccep.failrate <- lm(log10(ccep)~fail.rate+log10(D)+gm.stim+gm.rec,data=df.plt.pt)
          results$p[[wave]]$ccep.fail.rate[[pt]] <- get.coef.p.val(m.ccep.failrate,'fail.rate')
          results$b[[wave]]$ccep.fail.rate[[pt]] <- get.coef.beta(m.ccep.failrate,'fail.rate')
          if(identical(df.plt.pt$fail.rate,m.ccep.failrate$model$fail.rate)){
            df.plt.pt$ccep.fail.rate.pr <- get.partial.resids(m.ccep.failrate,'fail.rate')$y
          } else{error('Model and data frame are not aligned')} # put in a little check that everything is lining up
        } else{
          
          results$p[[wave]]$ccep.fail.rate[[pt]] <- results$b[[wave]]$ccep.fail.rate[[pt]] <- NA
          df.plt.pt$ccep.fail.rate.pr <- NA
          
        }
    } else{
      results$p[[wave]]$ccep.fail.rate[[pt]] <- results$b[[wave]]$ccep.fail.rate[[pt]] <- NA
      results$p[[wave]]$ccep.coefvar[[pt]] <- results$b[[wave]]$ccep.coefvar[[pt]] <- NA
      results$p[[wave]]$ccep.std[[pt]] <- results$b[[wave]]$ccep.std[[pt]] <- NA
    }
    df.plt.pt.all[[wave]][[pt]] <- df.plt.pt
  }
  
}

lapply(results$p,as.data.frame)
lapply(results$b,as.data.frame)

df.p.all <- list()
df.b.all <- list()
for(test in names(results$p$N1)){
  df.p <- as.data.frame(sapply(results$p, function(X) X[[test]]))
  df.p$pt <- rownames(df.p)
  df.p <- collapse.columns(df.p,cnames = waveforms,groupby='pt')
  df.p.all[[test]] <- rename.columns(df.p,c('names','group'),c('wave','pt'))
  
  df.b <- as.data.frame(sapply(results$b, function(X) X[[test]]))
  df.b$pt <- rownames(df.b)
  df.b <- collapse.columns(df.b,cnames = waveforms,groupby='pt')
  df.b.all[[test]] <- rename.columns(df.b,c('names','group'),c('wave','pt'))
  
}

# plot results of CCEP level analysis of spike rate vs. rho
# for each patient
plot.cfg <- list(list(ttl='CCEPVsCoefvar',
                      m.name='ccep.coefvar',
                      y.name='ccep.coefvar.pr',
                      x.name='coefvar'),
                 list(ttl='CCEPVsStd',
                      m.name='ccep.std',
                      y.name='ccep.std.pr',
                      x.name='std'),
                 list(ttl='CCEPVsFailRate',
                      m.name='ccep.fail.rate',
                      y.name='ccep.fail.rate.pr',
                      x.name='fail.rate'),
                 list(ttl='CCEPVsRho',
                      m.name='ccep.rho',
                      y.name='ccep.rho.pr',
                      x.name='rho'))

for(cfg in plot.cfg){
  for(wave in waveforms){
    
    df.plt.wave <- df.plt[df.plt$wave == wave,]
    df.plt.wave <- do.call(rbind,df.plt.pt.all[[wave]])
    df.p.wave <- df.p.all[[cfg$m.name]][df.p.all[[cfg$m.name]]$wave==wave,]
    
    p <- ggplot(df.plt.wave,aes_string(x=cfg$x.name,y=cfg$y.name)) + facet_wrap(~pt,scales = 'free')+ theme_minimal() + ggtitle(wave) +
      geom_point(alpha=0.5,stroke=0,size=0.5) + geom_smooth(method='lm',size=0.5) + standard_plot_addon() +
      geom_text(data=df.p.wave,aes(label=paste0('p = ',signif(values,2)),y=Inf,x=Inf),vjust=1,hjust=1,size=1)
    ggsave(plot = p,filename = paste0(savedir,cfg$ttl,wave,'.pdf'),width = 18,height=12,units= 'cm',useDingbats=FALSE)
    
  }
}

# plot results of CCEP level analysis of spike rate vs. rho
plot.cfg <- list(list(ttl='CCEPVsCoefvar',
                      m.name='ccep.coefvar',
                      y.name='ccep.coefvar.pr',
                      x.name='coefvar',
                      xl = 'Coefficient of Variation',
                      yl = 'log10(Amp.)'),
                 list(ttl='CCEPVsStd',
                      m.name='ccep.std',
                      y.name='ccep.std.pr',
                      x.name='std',
                      xl = 'Standard Deviation Across Trials',
                      yl = 'log10(Amp.)'),
                 list(ttl='CCEPVsFailRate',
                      m.name='ccep.fail.rate',
                      y.name='ccep.fail.rate.pr',
                      x.name='fail.rate',
                      xl = 'Failure Rate',
                      yl = 'log10(Amp.)'),
                 list(ttl='CCEPVsRho',
                      m.name='ccep.rho',
                      y.name='ccep.rho.pr',
                      x.name='rho',
                      xl = 'Rho',
                      yl = 'log10(Amp.)')
                 )
pt.plot <- 'HUP223'

for(cfg in plot.cfg){
  for(wave in waveforms){
    
    df.plt.wave <- df.plt.pt.all[[wave]][[pt.plot]]
    pval <- df.p.all[[cfg$m.name]][df.p.all[[cfg$m.name]]$wave==wave & df.p.all[[cfg$m.name]]$pt == pt.plot,'values']
    beta <- df.b.all[[cfg$m.name]][df.b.all[[cfg$m.name]]$wave==wave & df.b.all[[cfg$m.name]]$pt == pt.plot,'values']
    
    lab.p <- paste0(expression('p == '),scientific_10(signif(pval,2))) 
    lab.b <- paste0(expression(~beta),expression(' == '),scientific_10(signif(beta,2)))
    p <- ggplot(df.plt.wave,aes_string(x=cfg$x.name,y=cfg$y.name)) + theme_classic() + #ggtitle(wave) +
      geom_point(alpha=0.5,stroke=0,size=0.5,color=wes_palettes$BottleRocket2[2]) + 
      #stat_binhex() +
      #stat_density_2d(aes(fill=stat(density)), geom = 'raster', contour = FALSE) + scale_fill_viridis_c() +
      #coord_cartesian(expand = FALSE) +
      #nice_cbar('bottom') + 
      geom_smooth(method='lm',size=0.25,color=wes_palettes$BottleRocket2[3]) + 
      #geom_text(label=lab.p,y=Inf,x=Inf,vjust=1,hjust=1,size=1.5,parse=T) +
      #geom_text(label=lab.b,y=Inf,x=Inf,vjust=2,hjust=1,size=1.5,parse=T) +
      xlab(cfg$xl) + ylab(cfg$yl)+
      standard_plot_addon()
    ggsave(plot = p,filename = paste0(savedir,cfg$ttl,pt.plot,wave,'.pdf'),width = 4,height=4,units= 'cm',useDingbats=FALSE)
    
  }
}
 
