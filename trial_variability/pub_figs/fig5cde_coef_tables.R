# this script is same as Fig5d but controls for anatomic location, GM/WM status of electrodes and inter electrode distance 
# in assessing relationship between CCEP amplitude and SOZ

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
  
  df$Hippocampus.Stim <- df$BN.stim %in% c(215:218) | df$stim.name.ana == 'HIPP'
  df$Hippocampus.Rec <- df$BN.rec %in% c(215:218) | df$rec.name.ana == 'HIPP'
  
  df$Amygdala.Stim <- df$BN.stim %in% c(211:214) | df$stim.name.ana == 'AM'
  df$Amygdala.Rec <- df$BN.rec %in% c(211:214) | df$rec.name.ana == 'AM'
  
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
coef.names <- c('D','gm.stim','gm.rec','soz.non.sz')


# define functions to make regression table:

format.coefs <- function(x,sci=TRUE){
  x <- signif(x,2)
  sci.mask <- abs(log(abs(x),base=10))>2
  if(any(sci.mask)){
    x[sci.mask] <- sapply(x[sci.mask], function(y) format(y, scientific = sci))
  }
  return(x)
}

assemble.reg.tab <- function(m,dep.var,pt,wave){
  
  # nice names for model coefficient variables
  var.names <- c('(Intercept)','Distance','Gray Matter (Stim)','Gray Matter (Record)','Electrode In SOZ')
  
  # get coefficient table from lm summary object and format it in scientific/sig fig notation
  coef.tab <- summary(m)$coef
  coef.tab <- format.coefs(coef.tab)
  
  # add lines above for headers
  reg.tab <- rbind('','','',coef.tab,'','')
  
  # add clean headers for names of statistical parameters
  reg.tab[3,] <- c('Beta','Std. Error','t value','p value')
  
  # add f stat info
  f.stat.lines <- capture.output(print(summary(m)))
  reg.tab[nrow(reg.tab),1] <- f.stat.lines[grep('F-statistic',f.stat.lines)]
  
  # add clean variable names 
  reg.tab <- cbind('',reg.tab)
  #print(reg.tab)
  reg.tab[,1] <- c('','','',var.names,'','')
  
  # add headers
  reg.tab[1,1] <- paste0('Patient: ',pt)
  reg.tab[1,2] <- paste0('Wave: ',wave)
  reg.tab[1,3] <- paste0('Dependent Variable: ',dep.var)
  
  # remove row names
  rownames(reg.tab) <- NULL
  colnames(reg.tab) <- NULL
  return(reg.tab)
}

for(wave in waveforms){
  reg.tab.rho <- reg.tab.ccep <- reg.tab.std <- list()
  for(pt in pts){
    df.plt.pt <- df.plt[df.plt$pt==pt & !is.na(df.plt$soz.stim) & df.plt$wave == wave,]
    
    if(nrow(df.plt.pt)>0){
      if(length(unique(df.plt.pt$soz.non.sz[!is.na(df.plt.pt$ccep)]))>1){
        
        m.ccep.sozinout <- lm(log10(ccep)~log10(D)+gm.stim+gm.rec+soz.non.sz,data=df.plt.pt)

      } 
      
      m.rho.sozinout <- lm(rho~log10(D)+gm.stim+gm.rec+soz.non.sz,data=df.plt.pt)
      m.std.trials.sozinout <- lm(std.trials~log10(D)+gm.stim+gm.rec+soz.non.sz,data=df.plt.pt)
      m.log.std.trials.sozinout <- lm(log10(std.trials)~log10(D)+gm.stim+gm.rec+soz.non.sz,data=df.plt.pt)
      
      reg.tab.rho[[pt]] <- assemble.reg.tab(m.rho.sozinout,'Rho',pt,wave)
      reg.tab.rho[[pt]] <- cbind(reg.tab.rho[[pt]],'') # add blank row to horizontally concatenate
      reg.tab.ccep[[pt]] <- assemble.reg.tab(m.ccep.sozinout,'CCEP Amplitude',pt,wave)
      reg.tab.ccep[[pt]] <- cbind(reg.tab.ccep[[pt]],'') # add blank row to horizontally concatenate
      reg.tab.std[[pt]] <- assemble.reg.tab(m.std.trials.sozinout,'CCEP Amplitude Standard Deviation',pt,wave)
      reg.tab.std[[pt]] <- cbind(reg.tab.std[[pt]],'') # add blank row to horizontally concatenate
      
    }
  }
  
  reg.tab.rho <- do.call(cbind,reg.tab.rho)
  write.table(reg.tab.rho,file = paste0(savedir,'RegressionTable_',wave,'RhoSOZ_Fig5abc.csv'),row.names = F,col.names = F,sep=',')
  reg.tab.ccep <- do.call(cbind,reg.tab.ccep)
  write.table(reg.tab.ccep,file = paste0(savedir,'RegressionTable_',wave,'CCEPAmplitudeSOZ_Fig5abc.csv'),row.names = F,col.names = F,sep=',')
  reg.tab.std <- do.call(cbind,reg.tab.std)
  write.table(reg.tab.std,file = paste0(savedir,'RegressionTable_',wave,'CCEPStandardDeviationSOZ_Fig5abc.csv'),row.names = F,col.names = F,sep=',')
  
}
