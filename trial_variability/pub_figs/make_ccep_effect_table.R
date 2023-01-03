rm(list=ls())

basedir <- '/Users/Eli/Dropbox/CNTProjects/CCEPs_Projects/cceps_single_trials/'
source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))

pt.tbl <- read.csv(file = paste0(locations$results_folder,'PatientCharacteristicsCCEPs.csv'),check.names = F)

pt.tbl <- pt.tbl[,c('HUP ID','Seizure Onset Zone')]
rownames(pt.tbl) <- pt.tbl$`HUP ID`

###############################################################################
## Make table showing numbers of electrodes, CCEPs, and significant effects ###
###############################################################################

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

identical(rownames(pt.tbl),names(stim.chs))
pt.tbl[,'Sampling Rate (Hz)'] <- as.integer(unlist(fs))
pt.tbl[,'Stim Electrodes'] <- sapply(stim.chs,sum) # number of stimulation electrodes (before any artifact removal)
pt.tbl[,'Rec. Electrodes'] <- sapply(rec.chs,sum) # number of viable recording electrodes (before any artifact removal)
pt.tbl[,'Stim Trials'] <- as.integer(30*sapply(stim.chs,sum)) # Number of stimulation trials

# count non-artifactual CCEPs we measured - you can't analyze a CCEP where either N1 or N2 is 0 (b/c findpeaks couldn't find peak)
# but, the total number of nonartifactual CCEPs is greater than the CCEPs we analyzed (almost double) 
# Right now our code is written such that 0's come from either subthreshold CCEPS OR when findpeaks can't find a peak

pt.tbl[,'Suprathreshold CCEPs'] <- sapply(pt.tbl$`HUP ID`,function(pt) 
  sum(ave.cceps.N1[[pt]] > 0 & ave.cceps.N2[[pt]] > 0, na.rm = T))

# tile soz indicator to matrix size so you know if SOZ was stimulated for each CCEP
soz.all.pts <- unpack.mat.struct(elec.info,'soz')
soz.all.pts <- lapply(soz.all.pts, function(x) ifelse(x,yes='SOZ',no='Non-SOZ'))
soz.all.pts.rec <- lapply(soz.all.pts, function(x) matrix(x,nrow=length(x),ncol=length(x))) # check: colSums(soz.all.pts.rec$HUP211=='SOZ') should be 25 for every electrode
soz.all.pts.stim <- lapply(soz.all.pts, function(x) t(matrix(x,nrow=length(x),ncol=length(x)))) # check: rowSums(soz.all.pts.stim$HUP211=='SOZ') should be 25 for every electrode

waveforms <- c('N1','N2')
df.all <- list()
pct.sig.elecs.per.pt <- ct.sig.elecs.all.pts <- pct.sig.cceps.per.rec.elec.pt <- list()

for(wave in waveforms){
  print(wave)
  # unpack mat struct
  data.wave <- unpack.mat(spear.results,wave)
  spear.all.pts <- unpack.mat.struct(data.wave,'spear.trials')
  good.cceps.all.pts <- unpack.mat.struct(data.wave,'good.cceps')
  p.adj.all.pts <- unpack.mat.struct(data.wave,'spear.trials.p.adj')
  ave.cceps.wave <- unpack.mat.struct(ave.cceps,wave)
  
  chLabels.ana.all <- lapply(unpack.mat.struct(elec.info,'chLabels.ana'),cell.to.vec)
  
  # matrix whose elements correspond to name of stimulating electrode for each matrix element (i.e. CCEP)
  rec.name.ana <- lapply(chLabels.ana.all, function(x) matrix(x,nrow=length(x),ncol=length(x)))
  stim.name.ana <- lapply(chLabels.ana.all, function(x) t(matrix(x,nrow=length(x),ncol=length(x))))
  
  pt.indicator <- lapply(pts, function(pt) matrix(pt,nrow=nrow(spear.all.pts[[pt]]),ncol=ncol(spear.all.pts[[pt]])))

  df <- data.frame(rho=list.mat.to.vec(spear.all.pts),
                   g=list.mat.to.vec(good.cceps.all.pts),
                   p.adj=list.mat.to.vec(p.adj.all.pts),
                   stim.name.ana=list.mat.to.vec(stim.name.ana),
                   rec.name.ana=list.mat.to.vec(rec.name.ana),
                   pt=list.mat.to.vec(pt.indicator),
                   soz.stim=list.mat.to.vec(soz.all.pts.stim),
                   soz.rec=list.mat.to.vec(soz.all.pts.rec),
                   ccep=list.mat.to.vec(ave.cceps.wave),
                   wave=wave,
                   stringsAsFactors = F)
  
  df$soz.non.sz <- ifelse(df$soz.stim == 'Non-SOZ' & df$soz.rec == 'Non-SOZ',yes = 'Outside',no='Inside')
  df$soz.stim.rec <- paste0(df$soz.stim,' -> ',df$soz.rec)
  df$sig <- df$p.adj < 0.05
  
  # count number of significant effects at recording electrodes
  ct.sig.elecs.all.pts[[wave]] <- lapply(p.adj.all.pts,function(X) rowSums(X<0.05,na.rm=T))
  # 100 * number of recording electrodes with any significant monotonic trends / number of recording electrodes with any good cceps
  # i.e. how many recording electrodes with good cceps have this effect at all
  pct.sig.elecs.per.pt[[wave]] <- sapply(pts, function(pt) 100*sum(ct.sig.elecs.all.pts[[wave]][[pt]]>0) / sum(rowSums(good.cceps.all.pts[[pt]]) > 0))
  # 100 * number of significant effect at each recording electrode / number of CCEPs at each recording electrode
  # i.e. for each recording electrode how common is this effect
  pct.sig.cceps.per.rec.elec.pt[[wave]] <- sapply(pts, function(pt) 100*(ct.sig.elecs.all.pts[[wave]][[pt]] / rowSums(good.cceps.all.pts[[pt]])))
  
  df.good <- df[df$g==1,] # work with smaller df with no NAs to improve speed 
  df.all[[wave]] <- df
}

# Count number of GOOD (suprathreshold and NO artifact for both N1 and N2) CCEPs relative to the SOZ

good.ccep.mask <- df.all$N1$ccep >0 & df.all$N2$ccep >0
pt.tbl[,'SOZ -> SOZ'] <- sapply(pts, function(pt) sum(good.ccep.mask & df.all$N1$soz.stim == 'SOZ' & df.all$N1$soz.rec == 'SOZ' & df.all$N1$pt == pt,na.rm=T))
pt.tbl[,'SOZ -> Non-SOZ'] <- sapply(pts, function(pt) sum(good.ccep.mask & df.all$N1$soz.stim == 'SOZ' & df.all$N1$soz.rec == 'Non-SOZ' & df.all$N1$pt == pt,na.rm=T))
pt.tbl[,'Non-SOZ -> SOZ'] <- sapply(pts, function(pt) sum(good.ccep.mask & df.all$N1$soz.stim == 'Non-SOZ' & df.all$N1$soz.rec == 'SOZ' & df.all$N1$pt == pt,na.rm=T))
pt.tbl[,'Non-SOZ -> Non-SOZ'] <- sapply(pts, function(pt) sum(good.ccep.mask & df.all$N1$soz.stim == 'Non-SOZ' & df.all$N1$soz.rec == 'Non-SOZ' & df.all$N1$pt == pt,na.rm=T))

# Count the number of rejected (either for findpeaks or amplitude) CCEPs relative to SOZ
# good.ccep.mask <- df.all$N1$ccep ==0 | df.all$N2$ccep ==0
# pt.tbl[,'SOZ -> SOZ (Rej)'] <- sapply(pts, function(pt) sum(good.ccep.mask & df.all$N1$soz.stim == 'SOZ' & df.all$N1$soz.rec == 'SOZ' & df.all$N1$pt == pt,na.rm=T))
# pt.tbl[,'SOZ -> Non-SOZ (Rej)'] <- sapply(pts, function(pt) sum(good.ccep.mask & df.all$N1$soz.stim == 'SOZ' & df.all$N1$soz.rec == 'Non-SOZ' & df.all$N1$pt == pt,na.rm=T))
# pt.tbl[,'Non-SOZ -> SOZ (Rej)'] <- sapply(pts, function(pt) sum(good.ccep.mask & df.all$N1$soz.stim == 'Non-SOZ' & df.all$N1$soz.rec == 'SOZ' & df.all$N1$pt == pt,na.rm=T))
# pt.tbl[,'Non-SOZ -> Non-SOZ (Rej)'] <- sapply(pts, function(pt) sum(good.ccep.mask & df.all$N1$soz.stim == 'Non-SOZ' & df.all$N1$soz.rec == 'Non-SOZ' & df.all$N1$pt == pt,na.rm=T))

sig.eff.cts <- lapply(ct.sig.elecs.all.pts, function(X) sapply(X,sum))
sig.eff.pcts <- lapply(df.all, function(df)
  sapply(pts, function(pt) 100*mean(df[df$pt==pt & df$g==1,'sig'])))
lapply(sig.eff.pcts,summary)
pt.tbl[,'Monotonic Trends - N1 (%)'] <- paste0(sig.eff.cts$N1,' (',signif(sig.eff.pcts$N1,2),')')
pt.tbl[,'Monotonic Trends - N2 (%)'] <- paste0(sig.eff.cts$N2,' (',signif(sig.eff.pcts$N2,2),')')

# analyze phase analysis

data.effects <- readMat(paste0(locations$results_folder,'pub_figs/fig4/AllPatientsHippocampalPhaseCLcorr.mat'))
clcorr.results <- unpack.mat(data.effects,'clcorr.results')

for(wave in waveforms){
  data.wave <- unpack.mat(clcorr.results,wave)
  pts <- name(names(data.wave))
  pts.elecs <- lapply(pts, function(pt) names(unpack.mat(data.wave,pt)))
  clcorr.all.pts <- unpack.mat.struct.variable(data.wave,pts,pts.elecs,'clcorr.trials')
  p.adj.all.pts <- unpack.mat.struct.variable(data.wave,pts,pts.elecs,'clcorr.trials.p.adj')
  
  print(sapply(pts, function(pt) sum(do.call(rbind,p.adj.all.pts[[pt]])<0.05,na.rm=T)))
  df.plt <- data.frame(r=list.mat.to.vec(list.mat.to.vec(clcorr.all.pts)),
                       #g=list.mat.to.vec(good.cceps.all.pts),
                       p=list.mat.to.vec(list.mat.to.vec(p.adj.all.pts)))
  
  # count frequency of phase effect
  sig.count.all.pts <- lapply(p.adj.all.pts, function(X) sapply(X,function(X) which(X<0.05)))
  sig.count.all.pts <- lapply(sig.count.all.pts, function(X) length(unique(unlist(X))))
  good.cceps.all.pts <- lapply(p.adj.all.pts, function(X) sapply(X,function(X) which(!is.na(X))))
  good.cceps.all.pts <- lapply(good.cceps.all.pts, function(X) length(unlist(X)))
  print(paste(wave,'Hippocampal Phase: percent of CCEPs with significant effects'))
  sig.eff.pcts.all.pts <- sapply(pts, function(pt) 100*sig.count.all.pts[[pt]] / good.cceps.all.pts[[pt]])
  print(summary(sig.eff.pcts.all.pts))
  
  pt.tbl[,paste0('Hippocampal Phase - ',wave,' (%)')] <- paste0(unlist(sig.count.all.pts),' (',signif(sig.eff.pcts.all.pts,2),')')
  
}

print(xtable(pt.tbl,caption='CCEP measurement and prevalence of effects',label='table:tableS1'),include.rownames=FALSE)
