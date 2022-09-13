rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig3/fig3a/')
dir.create(savedir,recursive = T)

# load distribution of spearman correlations between trial index and N1/N2/pre stim negative control

# load distribution of spearman correlations between trial index and N1/N2/pre stim negative control

data <- readMat(paste0(locations$results_folder,'pub_figs/AllPatientsSpearmanLocsSoz.mat'))
spear.results <- unpack.mat(data,'spear.results')
elec.info <- unpack.mat(data,'elec.info')
pts <- name(names(elec.info))

# loop through wave forms, put all correlations and p values in one data frame to make one plot to show all distributions
# start by getting one data frame for each wave form (N1, N2, pre stim control) then combine at the end

waveforms <- names(spear.results)
df.all <- list()
pct.sig.elecs.per.pt <- list() # count number of RECORDING electrodes with significant effects per patient for each wave
for(wave in waveforms){
  
  # unpack mat struct
  data.wave <- unpack.mat(spear.results,wave)
  spear.all.pts <- unpack.mat.struct(data.wave,'spear.trials')
  good.cceps.all.pts <- unpack.mat.struct(data.wave,'good.cceps')
  p.adj.all.pts <- unpack.mat.struct(data.wave,'spear.trials.p.adj')
  D.all.pts <- unpack.mat.struct(elec.info,'D')
  pt.indicator <- lapply(pts, function(pt) matrix(pt,nrow=nrow(spear.all.pts[[pt]]),ncol=ncol(spear.all.pts[[pt]])))
  
  # get distribution of spearman correlations concatenated across all patients
  
  df.plt <- data.frame(r=list.mat.to.vec(spear.all.pts),
                       g=list.mat.to.vec(good.cceps.all.pts),
                       p=list.mat.to.vec(p.adj.all.pts),
                       d=list.mat.to.vec(D.all.pts),
                       pt=list.mat.to.vec(pt.indicator))
  
  # remove undetectable CCEPs
  df.plt <- df.plt[df.plt$g==1,]
  df.plt$wave <- wave
  
  # make indicator for adjusted p value < 0.05
  df.plt$sig <- df.plt$p < 0.05
  
  # plot histogram
  
  # p <- ggplot(df.plt) + geom_histogram(aes(x=r),alpha=0.5,fill = wes_palettes$BottleRocket1[4]) + theme_classic() +
  #   scale_y_continuous(expand = c(0,0)) + scale_x_continuous(limits=c(-1,1)) +
  #   #scale_fill_manual(values=wes_palettes$BottleRocket1) +
  #   xlab(expression("Spearman\'s"~rho)) + ylab('# of CCEPs') +
  #   theme(text=element_text(size=6))
  # count number of electrodes with significant effects
  
  ct.sig.elecs.all.pts <- lapply(p.adj.all.pts,function(X) rowSums(X<0.05,na.rm=T))
  pct.sig.elecs.per.pt[[wave]] <- sapply(pts, function(pt) 100*sum(ct.sig.elecs.all.pts[[pt]]>0) / sum(rowSums(good.cceps.all.pts[[pt]]) > 0))
  
  # store each data frame in a list 
  df.all[[wave]] <- df.plt
  
}

# compute the prevalence of this effect in each patient, for each wave
sig.eff.pcts <- lapply(df.all, function(df)
  sapply(pts, function(pt) 100*mean(df[df$pt==pt,'sig'])))
lapply(sig.eff.pcts,summary)

# combine the data frames into one list
df.all <- do.call(rbind,df.all)

## FIG 3A: Distribution of sequential effects in all patients -- boxplot for p > 0.05 to reduce overplot ##

p <- ggplot() + 
  ggdist::stat_halfeye(aes(x = df.all$wave, y = df.all$r),
                       adjust = .5, 
                       width = .3, 
                       .width = 0, 
                       justification = -.5, 
                       point_colour = NA,
                       alpha=0.8,
                       fill=wes_palettes$Royal1[4]
  ) + 
  geom_boxplot(aes(x=df.all$wave[!df.all$sig],y=df.all$r[!df.all$sig]),
               fill = wes_palettes$Royal1[1],
               alpha = 0.3,
               width = 0.1,
               size=0.1,
               outlier.size = 0.2,
               outlier.stroke = 0) +
  geom_point(aes(x=df.all$wave[df.all$sig],y=df.all$r[df.all$sig]),
             size = 0.2,
             alpha = 1,
             color = wes_palettes$Royal1[2],
             stroke=0,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) + 
  theme_bw() +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  ylab(expression("Spearman\'s"~rho)) + xlab('') + ggtitle('Monotonic Trends') +
  scale_y_continuous(limits=c(-1,1)) +
  standard_plot_addon() +
  theme(axis.title.y=element_text(size=6),legend.position='none',legend.box = 'none',legend.background = element_blank(),
        legend.margin = ggplot2::margin(0,0,0,0),
        legend.key.height = unit(0.01,units='cm'),legend.key.width=unit(1.25,'cm'))

ggsave(plot = p,filename = paste0(savedir,'Fig3a_left.pdf'),width = 4,height=4,units= 'cm',useDingbats=FALSE)

## FIG 3A: Distribution of sequential effects in all patients -- boxplot for p > 0.05 to reduce overplot ##

p <- ggplot() + 
  ggdist::stat_halfeye(aes(x = df.all$wave, y = df.all$d),
                       adjust = .5, 
                       width = .3, 
                       .width = 0, 
                       justification = -.5, 
                       point_colour = NA,
                       alpha=0.8,
                       fill=wes_palettes$Royal1[4]
  ) + 
  geom_boxplot(aes(x=df.all$wave[!df.all$sig],y=df.all$d[!df.all$sig]),
               fill = wes_palettes$Royal1[1],
               alpha = 0.3,
               width = 0.1,
               size=0.1,
               outlier.size = 0.2,
               outlier.stroke = 0,
               outlier.alpha = 0.2) +
  geom_point(aes(x=df.all$wave[df.all$sig],y=df.all$d[df.all$sig]),
             size = 0.2,
             alpha = 1,
             color = wes_palettes$Royal1[2],
             stroke=0,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) + 
  theme_bw() +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  ylab("Inter-Electrode Distance (mm)") + xlab('') + ggtitle('Spatial Distribution') +
  standard_plot_addon() +
  theme(axis.title.y=element_text(size=6),legend.position='none',legend.box = 'none',legend.background = element_blank(),
        legend.margin = ggplot2::margin(0,0,0,0),
        legend.key.height = unit(0.01,units='cm'),legend.key.width=unit(1.25,'cm'))

ggsave(plot = p,filename = paste0(savedir,'Fig3a_right.pdf'),width = 4,height=4,units= 'cm',useDingbats=FALSE)

## Code below has full jitter for the distribution which makes an unwieldy vector graphic

## FIG 3A: Distribution of sequential effects in all patients ##

# p <- ggplot(df.all, aes(x = wave, y = r)) + 
#   ggdist::stat_halfeye(
#     adjust = .5, 
#     width = .3, 
#     .width = 0, 
#     justification = -.5, 
#     point_colour = NA,
#     alpha=0.8,
#     fill=wes_palettes$Royal1[4]
#   ) + 
#   geom_point(aes(color=sig,alpha=sig),
#              size = 0.3,
#              #alpha = .8,
#              stroke=0,
#              position = position_jitter(
#                seed = 1, width = .1
#              )
#   ) + 
#   scale_color_manual(values=wes_palettes$Royal1[c(1,2)],breaks = c(FALSE,TRUE),labels= c(TeX('$p_{FDR}>0.05$'),TeX('$p_{FDR}$<0.05')),name='') +
#   scale_alpha_manual(values=c(0.1,1),breaks=c(FALSE,TRUE),limits=c(FALSE,TRUE),guide='none')+
#   theme_bw() +
#   coord_cartesian(xlim = c(1.2, NA), clip = "off") +
#   ylab(expression("Spearman\'s"~rho)) + xlab('') + 
#   standard_plot_addon() +
#   theme(axis.title.y=element_text(size=6),legend.position='none',legend.box = 'none',legend.background = element_blank(),
#         legend.margin = ggplot2::margin(0,0,0,0),
#         legend.key.height = unit(0.01,units='cm'),legend.key.width=unit(1.25,'cm'))
# 
# ggsave(plot = p,filename = paste0(savedir,'Fig3a.pdf'),width = 4,height=4,units= 'cm')
# 
# ## FIG 3B: Inter-electrode Distance for CCEPs stratified by significance of sequential effects ##
# 
# p <- ggplot(df.all, aes(x = wave, y = d)) + 
#   ggdist::stat_halfeye(
#     adjust = .5, 
#     width = .3, 
#     .width = 0, 
#     justification = -.5, 
#     point_colour = NA,
#     alpha=0.8,
#     fill=wes_palettes$Royal1[4]
#   ) + 
#   geom_point(aes(color=sig,alpha=sig),
#              size = 0.3,
#              #alpha = .8,
#              stroke=0,
#              position = position_jitter(
#                seed = 1, width = .1
#              )
#   ) + 
#   scale_color_manual(values=wes_palettes$Royal1[c(1,2)],breaks = c(FALSE,TRUE),labels= c(TeX('$p_{FDR}>0.05$'),TeX('$p_{FDR}$<0.05')),name='') +
#   scale_alpha_manual(values=c(0.1,1),breaks=c(FALSE,TRUE),limits=c(FALSE,TRUE),guide='none')+
#   theme_bw() +
#   coord_cartesian(xlim = c(1.2, NA), clip = "off") +
#   ylab('Inter-Electrode Distance (mm)') + xlab('') + 
#   standard_plot_addon() +
#   theme(axis.title.y =element_text(size=6),legend.position='none',legend.box = 'none',legend.background = element_blank(),
#         legend.margin = ggplot2::margin(0,0,0,0),
#         legend.key.height = unit(0.01,units='cm'),legend.key.width=unit(1.25,'cm'))
# 
# ggsave(plot = p,filename = paste0(savedir,'Fig3b.pdf'),width = 4,height=4,units= 'cm')
# 
