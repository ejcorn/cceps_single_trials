# this script is same as Fig5d but controls for anatomic location, GM/WM status of electrodes and inter electrode distance 
# in assessing relationship between CCEP amplitude and SOZ

rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/miscfxns/statfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig5/')
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

chLabels.all <- lapply(unpack.mat.struct(elec.info,'chLabels'),cell.to.vec)
pt.indicator.elecs <- lapply(pts, function(pt) rep(pt,length(chLabels.all[[pt]])))

# plot spikes vs. SOZ

df.ss <- data.frame(soz = list.mat.to.vec(soz.all.pts),
                    spikes = list.mat.to.vec(spikes.all.pts),
                    pt = list.mat.to.vec(pt.indicator.elecs))
df.ss <- df.ss[!is.na(df.ss$spikes) & !is.na(df.ss$soz),]
df.ss.pt <- lapply(pts, function(pt) df.ss[df.ss$pt==pt,])
p.vals <- sapply(df.ss.pt, function(df) if(nrow(df)>0){wilcox.test(formula=spikes~soz,data=df)$p.value} else{NA})
df.p <- data.frame(pt=pts,p=p.vals)
df.p <- df.p[!is.na(df.p$p),]
# paste0(expression('p ==\n'),scientific_10(signif(p,2))))

p <- ggplot(df.ss) + geom_boxplot(aes(x=pt,y=spikes,fill=soz),size=0.1,outlier.size=0.2,outlier.alpha = 0.5,outlier.stroke=0) + 
  theme_classic() + geom_text(data=df.p,aes(x=pt,y=Inf,label=p.signif(p)),vjust=1,size=1.5) +
  scale_fill_manual(values=wes_palettes$BottleRocket2[c(3,2)])+
  ylab('Spikes/min')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5)) +standard_plot_addon() +
  theme(legend.position = 'bottom',legend.key.size = unit(0.1,'cm'),
        legend.title = element_blank(),
        axis.title.x = element_blank())
ggsave(plot = p,filename = paste0(savedir,'Fig5a.pdf'),width = 4,height=4,units= 'cm',useDingbats=FALSE)
