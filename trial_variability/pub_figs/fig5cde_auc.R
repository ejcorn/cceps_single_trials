# this script uses caret to check predictive power of 

rm(list=setdiff(ls(),'basedir'))

source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))
source(paste0(basedir,'trial_variability/miscfxns/statfxns.R'))
source(paste0(basedir,'trial_variability/plottingfxns/plottingfxns.R'))

# set save directory
savedir <- paste0(locations$results_folder,'pub_figs/fig5/fig5cde_auc/')
dir.create(savedir,recursive = T)

load(paste0(savedir,'TrainData.RData'))

wave <- 'N1'
waveforms <- c('N1','N2')

createMultiGroupKFold <- function(group,k,times){
  # only works if times is <1000
  groupFolds <- lapply(1:times,function(x) groupKFold(group = group,k = k))
  for(ti in 1:times){
    names(groupFolds[[ti]]) <- paste0(names(groupFolds[[ti]]),'.Rep',sprintf('%03d',ti))
  }
  return(Reduce('c',groupFolds))
}

wave <- 'N1'
df.train <- df.all[[wave]]
groupFolds <- createMultiGroupKFold(df.train$pt,5,3)

ctrl <- trainControl(method = "repeatedcv",classProbs = T,savePredictions = T,
                     index=groupFolds,summaryFunction = twoClassSummary)

covariates <- 'log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+log10(ccep)+log10(std.trials)'
glmFit1 <- train(reformulate(covariates,response = 'soz.non.sz'), data = df.train, 
                    method = "glm", 
                    trControl = ctrl,
                  na.action = na.pass)

f <- formula(soz.non.sz.binary ~ log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec+log10(ccep)+log10(std.trials))
#f <- formula(soz.non.sz.binary ~ log10(D)+gm.stim+gm.rec+Somatosensory.Stim+Somatosensory.Rec+Motor.Stim+Motor.Rec+Visual.Stim+Visual.Rec)
m <- glm(f,data=df.train)
m.r <- roc(response=df.train$soz.non.sz.binary,predictor = predict(m,df.train,type='response'))
m.r$auc

df.pt <- df.train[df.train$pt == 'HUP223',]
m <- glm(f,data=df.pt)
m.r <- roc(response=df.pt$soz.non.sz.binary,predictor = predict(m,df.pt,type='response'))
m.r$auc

ctrl <- trainControl(method = "repeatedcv",classProbs = T,savePredictions = T,
                    summaryFunction = twoClassSummary)
glmFit.pt <- train(reformulate(covariates,response = 'soz.non.sz'), data = df.pt, 
                 method = "glm", 
                 trControl = ctrl,
                 na.action = na.pass)
ggplot(glmFit.pt$pred) + geom_histogram(aes(x=Inside,fill=obs),alpha=0.1,position = 'identity')

xgb_grid_1 <- expand.grid(
  nrounds = 10,
  eta = 0.3,
  max_depth = 5,
  gamma = 0,
  colsample_bytree=0.75, 
  min_child_weight=1,
  subsample=0.75
)

ctrl <- trainControl(method = "repeatedcv",number = 5,repeats = 3,classProbs = T,savePredictions = T)
xgbFit.pt <- train(reformulate(covariates,response = 'soz.non.sz'), data = df.pt, 
               method = "xgbTree", 
               trControl = ctrl,
               tuneGrid = xgb_grid_1,
               verbose = FALSE,
               importance = FALSE,
               metric='ROC',
               preProcess = c('center','scale'),
               nthread=1)

ggplot(xgbFit.pt$pred) + geom_histogram(aes(x=Inside,fill=obs),alpha=0.1,position = 'identity')

############


for(wave in waveforms){
  for(pt in pts){
    df.plt.pt <- df.plt[df.plt$pt==pt & !is.na(df.plt$soz.stim) & df.plt$wave == wave,]
      if(sum(df.plt.pt$soz.non.sz.binary)>0){
        m.ccep.sozinout <- lm(log10(ccep)~log10(D)+gm.stim+gm.rec+soz.non.sz,data=df.plt.pt)
        m.rho.sozinout <- lm(rho~log10(D)+gm.stim+gm.rec+soz.non.sz,data=df.plt.pt)
        m.std.trials.sozinout <- lm(std.trials~log10(D)+gm.stim+gm.rec+soz.non.sz,data=df.plt.pt)
        m.log.std.trials.sozinout <- lm(log10(std.trials)~log10(D)+gm.stim+gm.rec+soz.non.sz,data=df.plt.pt)
        
        f.base <- reformulate(termlabels = paste0(covariates,'') ,response='soz.non.sz.binary')
        m.auc$base <- boot.auc.df(df.plt.pt,f.base)
        f.all <- reformulate(termlabels = paste0(covariates,'+rho+log10(ccep)+std.trials') ,response='soz.non.sz.binary')
        m.auc$all <- boot.auc.df(df.plt.pt,f.all)
        f.ccep <- reformulate(termlabels = paste0(covariates,'+log10(ccep)') ,response='soz.non.sz.binary',intercept = T)
        m.auc$ccep <- boot.auc.df(df.plt.pt,f.ccep)
        f.std <- reformulate(termlabels = paste0(covariates,'+std.trials') ,response='soz.non.sz.binary')
        m.auc$std <- boot.auc.df(df.plt.pt,f.std)
        f.rho <- reformulate(termlabels = paste0(covariates,'+rho') ,response='soz.non.sz.binary')
        m.auc$rho <- boot.auc.df(df.plt.pt,f.rho)
    
        
        auc.cis <- as.data.frame(t(sapply(m.auc,function(X) X)))
        colnames(auc.cis) <- c('ymin','y','ymax')
        auc.cis$names <- rownames(auc.cis)
        

      }
  }
    df.plt.pt.all[[wave]][[pt]] <- df.plt.pt
    p.all <- plot_grid(plotlist=p.all[[wave]],
                       align='hv',ncol=3)
    ggsave(plot = p.all,filename = paste0(savedir,'Fig5cdeAUC',wave,'.pdf'),width = 18,height=12,units= 'cm',useDingbats=FALSE)
    
}


x.order <- c('base','all','ccep','rho','std')
p.all[[wave]][[pt]] <- ggplot(auc.cis) + geom_col(aes(x=names,y=y)) +
  geom_errorbar(aes(x=names,y=y,ymin=ymin,ymax=ymax)) +
  scale_x_discrete(limits=x.order,name='') + 
  scale_y_continuous(limits=c(0,1),name='AUC')+
  ggtitle(paste0(wave,': ',pt))+
  standard_plot_addon() + x_text_90()

lapply(results.soz,as.data.frame)
lapply(results.soz.perm,as.data.frame)
lapply(m.soz.pred.auc$N1,function(X) auc(X))

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
                 list(ttl='FailRateVsSOZ',
                      m.name='fail.rate.soz.inout',
                      y.name='fail.rate.soz.pr',
                      x.name='soz.non.sz',
                      xl='SOZ',
                      yl='Fail Rate'),
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
ggsave(plot = p.all,filename = paste0(savedir,'Fig5cde.pdf'),width = 8.75,height=10,units= 'cm',useDingbats=FALSE)
