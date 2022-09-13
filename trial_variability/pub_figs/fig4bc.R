rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig4/')
dir.create(savedir,recursive = T)

# load distribution of clcorr between trial index and N1/N2/pre stim negative control

data <- readMat(paste0(locations$results_folder,'pub_figs/AllPatientsSpearmanLocsSoz.mat'))
plot.list <- list(
  list(data.roi = readMat(paste0(locations$results_folder,'pub_figs/fig4/AllPatientsHippocampalPhaseCLcorr.mat')),
                          analysis.name = 'Hippocampal Phase',
                          fig.name = 'Fig4b',
                          metric='clcorr'),
  list(data.roi = readMat(paste0(locations$results_folder,'pub_figs/fig4/AllPatientsRecordingElectrodePhaseCLcorr.mat')),
       analysis.name = 'Recording Electrode Phase',
       fig.name = 'Fig4c',
       metric='clcorr'
  ),
     list(data.roi = readMat(paste0(locations$results_folder,'pub_figs/fig4/AllPatientsACCPhaseCLcorr.mat')),
                             analysis.name = 'Anterior Cingulate Phase',
                             fig.name = 'FigXX',
                             metric='clcorr'
          )#,
  # list(data.roi = readMat(paste0(locations$results_folder,'pub_figs/fig4/AllPatientsRecordingElectrodePowerSpearman.mat')),
  #      analysis.name = 'Recording Electrode Power',
  #      fig.name = 'Fig5c',
  #      metric='spear'
  # ),
  # list(data.roi = readMat(paste0(locations$results_folder,'pub_figs/fig4/AllPatientsHippocampalPowerSpearman.mat')),
  #      analysis.name = 'Hippocampal Power',
  #      fig.name = 'FigS_HippocampalPower',
  #      metric='spear'
  # )
)
 
for(pl in plot.list){
  print(pl$analysis.name)
  data.effects <- pl$data.roi
  clcorr.results <- unpack.mat(data.effects,paste0(pl$metric,'.results'))
  
  # loop through wave forms, put all correlations and p values in one data frame to make one plot to show all distributions
  # start by getting one data frame for each wave form (N1, N2, pre stim control) then combine at the end
  
  waveforms <- names(clcorr.results)
  df.all <- list()
  wave <- 'N2'
  for(wave in waveforms){
    
    # unpack mat struct
    data.wave <- unpack.mat(clcorr.results,wave)
    
    pts <- name(names(data.wave))
    pts.elecs <- lapply(pts, function(pt) names(unpack.mat(data.wave,pt)))
    
    if(pl$analysis.name == 'Recording Electrode Phase'){
      clcorr.all.pts <- unpack.mat.struct(data.wave,paste0(pl$metric,'.trials'))
      p.adj.all.pts <- unpack.mat.struct(data.wave,paste0(pl$metric,'.trials.p.adj'))
      #print(wave)
      #print(sapply(p.adj.all.pts,function(X) sum(X<0.05,na.rm=T)))
      df.plt <- data.frame(r=list.mat.to.vec(clcorr.all.pts),
                           #g=list.mat.to.vec(good.cceps.all.pts),
                           p=list.mat.to.vec(p.adj.all.pts))
    } else{
      clcorr.all.pts <- unpack.mat.struct.variable(data.wave,pts,pts.elecs,paste0(pl$metric,'.trials'))
      p.adj.all.pts <- unpack.mat.struct.variable(data.wave,pts,pts.elecs,paste0(pl$metric,'.trials.p.adj'))
      print(wave)
      print(sapply(pts, function(pt) sum(do.call(rbind,p.adj.all.pts[[pt]])<0.05,na.rm=T)))
      df.plt <- data.frame(r=list.mat.to.vec(list.mat.to.vec(clcorr.all.pts)),
                           #g=list.mat.to.vec(good.cceps.all.pts),
                           p=list.mat.to.vec(list.mat.to.vec(p.adj.all.pts)))
      
      # count frequency of phase effect
      sig.count.all.pts <- lapply(p.adj.all.pts, function(X) sapply(X,function(X) which(X<0.05)))
      sig.count.all.pts <- lapply(sig.count.all.pts, function(X) length(unique(unlist(X))))
      good.cceps.all.pts <- lapply(p.adj.all.pts, function(X) sapply(X,function(X) which(is.na(X))))
      good.cceps.all.pts <- lapply(good.cceps.all.pts, function(X) length(unique(unlist(X))))
      print(paste(wave,pl$analysis.name,': percent of CCEPs with significant effects'))
      print(summary(sapply(pts, function(pt) 100*sig.count.all.pts[[pt]] / good.cceps.all.pts[[pt]])))
      
      
    }
    # get distribution of phase-waveform correlations concatenated across all patients
    

    
    # remove undetectable CCEPs
    df.plt <- df.plt[!is.na(df.plt$r),]
    df.plt$wave <- wave
    
    # make indicator for adjusted p value < 0.05
    df.plt$sig <- df.plt$p < 0.05
    
    # plot histogram
    
    # p <- ggplot(df.plt) + geom_histogram(aes(x=r),alpha=0.5,fill = wes_palettes$BottleRocket1[4]) + theme_classic() +
    #   scale_y_continuous(expand = c(0,0)) + scale_x_continuous(limits=c(-1,1)) +
    #   #scale_fill_manual(values=wes_palettes$BottleRocket1) +
    #   xlab(expression("clcorrman\'s"~rho)) + ylab('# of CCEPs') +
    #   theme(text=element_text(size=6))
    # 
    
    # store each data frame in a list 
    
    df.all[[wave]] <- df.plt
    
  }
  
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
               size = 0.4,
               alpha = 1,
               color = wes_palettes$Royal1[2],
               stroke=0,
               position = position_jitter(
                 seed = 1, width = .1
               )
    ) + 
    theme_bw() +
    coord_cartesian(xlim = c(1.2, NA), clip = "off") +
    ylab(expression("Circular-Linear"~rho)) + xlab('') + ggtitle(pl$analysis.name)+
    standard_plot_addon() +
    theme(axis.title.y=element_text(size=6),
          legend.position='bottom',legend.box = 'none',
          #legend.background = element_blank(),
          legend.margin = ggplot2::margin(0,0,0,0),
          legend.key.height = unit(0.01,units='cm'),legend.key.width=unit(1.25,'cm'))
  
  p
  ggsave(plot = p,filename = paste0(savedir,pl$fig.name,'.pdf'),width = 4,height=4,units= 'cm',useDingbats=FALSE)
  
}  


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
  #   ylab(expression("clcorrman\'s"~rho)) + xlab('') + 
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
