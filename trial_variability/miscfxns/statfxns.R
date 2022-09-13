################
### p-values ###
################

p.signif <- function(p,ns='ns'){
  # take p values and make asterisks
  #   ns: p > 0.05
  # *: p <= 0.05
  # **: p <= 0.01
  # ***: p <= 0.001
  # ****: p <= 0.000001
  p.new <- rep('',length(p))
  p.new[is.na(p)] <- ''
  p.new[p > 0.05] <- ns
  p.new[p < 0.05 & p > 0.01] <- '*'
  p.new[p <= 0.01 & p > 0.001] <- '**'
  p.new[p <= 0.001 & p > 0.000001] <- '***'
  p.new[p <= 0.000001] <- '****'
  names(p.new) <- names(p)
  return(p.new)
  
}

p.signif.matrix <- function(p){
  # take matrix of p values and make asterisks
  #   ns: p > 0.05
  # *: p <= 0.05
  # **: p <= 0.01
  # ***: p <= 0.001
  # ****: p <= 0.000001
  p.new <- matrix(data = '',nrow = nrow(p),ncol=ncol(p))
  p.new[p > 0.05] <- 'ns'
  p.new[p < 0.05 & p > 0.01] <- '*'
  p.new[p <= 0.01 & p > 0.001] <- '**'
  p.new[p <= 0.001 & p > 0.000001] <- '***'
  p.new[p <= 0.000001] <- '****'
  return(p.new)
  
}

matrix.fdr.correct <- function(pvals){
 pvals.mat <- matrix(p.adjust(as.vector(as.matrix(pvals)), method='fdr'),ncol=ncol(pvals))
 colnames(pvals.mat) <- colnames(pvals)
 return(pvals.mat) 
}

list.vec.fdr.correct <- function(X){
  # for a list where each element is a vector of p-values
  # FDR correct over all p-values and return to list
  # correct
  p.fdr <- p.adjust(unlist(X),method='fdr')
  # get lengths --> cumulative indices
  X.ind <- cumsum(sapply(X,length))
  # initialize
  X.fdr <- list()
  X.names <- names(X)
  # start off list with first portion
  X.fdr[[X.names[1]]] <- p.fdr[1:X.ind[1]]
  # fill in list elements 2 to n
  for(i in 2:length(X.ind)){
    X.fdr[[X.names[i]]] <- p.fdr[(1+X.ind[i-1]):(X.ind[i])]
  }
  return(X.fdr)
}

list.fdr.correct <- function(X){
  # unlist a list, fdr correct over all values
  # relist the list in the same structure and return
  return(relist(flesh=p.adjust(unlist(X),method='fdr'),skeleton=X))
}

perm.test.matrix.2tail <- function(test.mat,dist.mat){
  
  # INPUTS:
  # test.mat: NxK matrix of test statistics
  # dist.mat: NxKxnperms matrix of null distribution of test statistics corresponding to test.mat
  #
  # OUTPUTS:
  # pvals.2tail: return 2 tailed p-values
  
  dims <- dim(dist.mat)
  nperms <- dims[3]
  test.mat <- abind(lapply(1:nperms, function(x) test.mat),along=3) # duplicate here to make array dimensions align
  pvals.2tail <- 2*pmin(apply(dist.mat >= test.mat,c(1,2),mean),apply(dist.mat <= test.mat,c(1,2),mean))
  return(pvals.2tail)
}


wilcox.auc <- function(pred,obs){
  # INPUTS:
  # pred: 
  # obs:
  #
  # OUTPUTS:
  #

  classes <- unique(obs)
  wt <- wilcox.test(pred[obs == classes[1]],pred[obs == classes[2]])
  return(wt$p.value)
}

get.partial.resids <- function(m,xname){
  # INPUTS:
  # m: lm object
  # xname: character of x variable name
  # xname.coefs: for factors when lm() appends level to coefficient name
  #
  # OUTPUTS:
  # df: with columns x and y to make partial residual plot
  # partial residual wrt x_i is (y-y^) + b_i*x_i
  
  pr <- m$residuals + m$model[,xname]*m$coefficients[xname] + m$coefficients['(Intercept)']
  df <- data.frame(x = m$model[,xname],y=pr) 
  return(df)
  
}

get.partial.resids <- function(m,xname){
  # INPUTS:
  # m: lm object
  # xname: character of x variable name
  # xname.coefs: for factors when lm() appends level to coefficient name
  #
  # OUTPUTS:
  # df: with columns x and y to make partial residual plot
  # partial residual wrt x_i is (y-y^) + b_i*x_i
  
  pr <- m$residuals + m$model[,xname]*m$coefficients[xname] + m$coefficients['(Intercept)']
  df <- data.frame(x = m$model[,xname],y=pr) 
  return(df)
  
}
get.partial.resids.factor <- function(m,xname,xnamelevel){
  # INPUTS:
  # m: lm object
  # xname: character of x variable name
  # xname.coefs: for factors when lm() appends level to coefficient name
  #
  # OUTPUTS:
  # df: with columns x and y to make partial residual plot
  # partial residual wrt x_i is (y-y^) + b_i*x_i
  
  pr <- m$residuals + as.numeric(m$model[,xname] == xnamelevel)*m$coefficients[paste0(xname,xnamelevel)] + m$coefficients['(Intercept)']
  df <- data.frame(x = m$model[,xname],y=pr)
  return(df)
  
}
#########################################################
### Permutation testing for CCEP electrode categories ###
#########################################################

perm.elec.cats <- function(cat.all.pts,nperms=1000){
    # INPUTS:
    # cat.all.pts: list with one element per patient of character vectors corresponding to electrode categories
    # OUTPUTS:
    # cat.all.pts.shuffled: cat.all.pts but with those character vectors shuffled randomly without replacement
    return(lapply(1:nperms,function(X) lapply(cat.all.pts, function(X) sample(X,replace=F))))
  }

perm.group.elec.test.stat <- function(metric.all.pts,cat.all.pts,fxn){
  # INPUTS:
  # metric.all.pts: list with one element per patient of numeric vectors corresponding to electrode metrics (i.e. in degree, number of significant spearmans, etc.)
  # cat.all.pts: list with one element per patient of character vectors corresponding to electrode categories
  # fxn: fxn applied to data frame with columns "metric" and "cat" for category to return test statistic
  # OUTPUTS:
  # cat.all.pts.shuffled: cat.all.pts but with those character vectors shuffled randomly without replacement
  df <- data.frame(metric=list.mat.to.vec(metric.all.pts),cat=list.mat.to.vec(cat.all.pts))
  return(fxn(df))
}

perm.test.elecs <- function(cat.all.pts,metric.all.pts,fxn,nperms=1000){
  # put it all together
  cat.all.pts.shuffle <- perm.elec.cats(cat.all.pts,nperms)  
  perm <- sapply(cat.all.pts.shuffle, function(X) perm.group.elec.test.stat(metric.all.pts,X,fxn))
  actual <- perm.group.elec.test.stat(metric.all.pts,cat.all.pts,fxn)
  results <- list(p=pval.2tail.np(0,perm-actual),
                  observed=actual,
                  null=perm)
  return(results)
}


lm.perm.null <- function(m,df.m,nperms=10){
  
  lm.form <- formula(m)
  response <- all.vars(lm.form)[1]
  df.m <- df.m[!is.na(df.m[,response]),]
  df.response.null <- lapply(1:nperms, function(x) df.m[sample(1:nrow(df.m)),response,drop=F])
  df.m.terms <- df.m[,-which(colnames(df.m) == response)]
  m.coefs <- coef(m)
  m.perm.diff <- lapply(df.response.null, function(x) coefficients(lm(formula=lm.form,data=cbind(x,df.m.terms))) - m.coefs)
  coefs.null <- do.call(rbind,m.perm.diff)
  pvals.all <- sapply(colnames(coefs.null), function(x) pval.2tail.np(0,coefs.null[,x]))
  return(pvals.all)
}

#########################################################################
### Permutation testing for CCEP electrode categories (Network Level) ###
#########################################################################

lm.perm.null.df.list <- function(m,df.null.list,coef.names){
  
  lm.form <- formula(m)
  m.coefs <- coef(m)
  m.perm.diff <- list()
  for(p in 1:length(df.null.list)){
    df.null <- df.null.list[[p]]
    if(all(sapply(df.null[complete.cases(df.null[,coef.names]),coef.names],function(x) length(unique(x))>1))){
      m.perm.diff[[p]] <- coefficients(lm(formula=lm.form,data=df.null)) - m.coefs
    } else{print(paste0(p,unique(df.null$pt)))}
  }
  #m.perm.diff <- lapply(df.null.list, function(df.null) coefficients(lm(formula=lm.form,data=df.null)) - m.coefs)
  coefs.null <- do.call(rbind,m.perm.diff)
  pvals.all <- sapply(colnames(coefs.null), function(x) pval.2tail.np(0,coefs.null[,x]))
  return(pvals.all)
}

#################
### Bootstrap ###
#################

meanboot <- function(data, i){
  d <- data[i]
  return(mean(d))   
}

medianboot <- function(data, i){
  d <- data[i]
  return(median(d))   
}

############################
### Mixed effects models ###
############################


lme.ms <- function(cfg,df){
  # INPUTS:
  # cfg: list with the named elements below
  #   fixed, random: character names of independent variables in fixed and random effects
  #   response: character name of response (default to 'Score' for this script)
  #   id: character name of ids (Default to 'scanid' for this script)
  # df: data frame containing fixed.IVs and random.IVs and response
  # OUTPUTS:
  # lme call with input fixed and random effects
  # use this so I can do easy model selection without copy-pasting the same parameters
  fixed <- reformulate(termlabels = cfg$fixed,response = cfg$response,intercept = T) # make fixed formula
  random <- reformulate(termlabels = cfg$random,response = cfg$response,intercept = T) # make random formula
  random <- as.formula(paste(Reduce(paste,deparse(random)),'|',cfg$id)) # add id variable to random
  ctrl <- lmeControl(opt = 'optim') 
  # parse formulas as text in order to extract coefficients and compare models with anova
  return(        eval(bquote(   lme(fixed=.(fixed),random=.(random), 
                                    data = df,na.action=na.exclude,control=ctrl,method='ML')   )) )
}

lme4.ms <- function(cfg,df){
  # INPUTS:
  # cfg: list with the named elements below
  #   fixed, random: character names of independent variables in fixed and random effects
  #   response: character name of response (default to 'Score' for this script)
  #   id: character name of ids (Default to 'scanid' for this script)
  # df: data frame containing fixed.IVs and random.IVs and response
  # OUTPUTS:
  # lme call with input fixed and random effects
  # use this so I can do easy model selection without copy-pasting the same parameters
  fixed <- reformulate(termlabels = cfg$fixed,response = cfg$response,intercept = T) # make fixed formula
  f <- as.formula(paste(cfg$response,'~',as.character(fixed)[3],'+(',paste(cfg$random,collapse='+'),'|',cfg$id,')'))
  ctrl <- lmerControl(optCtrl = list(method='Nelder-Mead')) 
  # parse formulas as text in order to extract coefficients and compare models with anova
  return(        eval(bquote(   lmer(formula = .(f), 
                                    data = df,na.action=na.exclude,control=ctrl,REML=F)   )) )
}

lme.compare <- function(mdl.gold,mdl.test,coef.sig){
  # INPUT:
  # mdl.gold, mdl.test: compare new model (mdl.test) to current gold std model (mdl.gold)
  # using 3 criteria that all must be true:
  # (c1) p < 0.05 for anova(m1,m2) AND 
  # (c2) AIC drops AND 
  # (c3) added coefficient is significant
  # coef.sig: name of coefficient to test significance for
  #
  # OUTPUT:
  # cond: true/false is model better
  c1 <- sort(anova(mdl.gold,mdl.test)$`p-value`) < 0.05 # likelihood ratio p < 0.05
  c2 <- summary(mdl.gold)$AIC > summary(mdl.test)$AIC # AIC decreases
  # print('1. Anova')
  # print(anova(mdl.gold,mdl.test))
  # print('2. AIC: gold vs test')
  # print(c(summary(mdl.gold)$AIC,summary(mdl.test)$AIC))
  print(paste('3. checking p-val of',coef.sig))
  print(summary(mdl.test)$tTable[coef.sig,,drop=F])
  c3 <- summary(mdl.test)$tTable[coef.sig,'p-value'] < 0.05 # p-value of most recently added coefficient is < 0.05
  # if the two models only differ by random effects, don't assess coefficient significance
  if(identical(names(mdl.gold$coefficients$fixed),names(mdl.test$coefficients$fixed))){c3 <- TRUE}
  cond <- all(c(c1,c2,c3)) # all conditions must be true to advance
  return(cond)
}

lme.last.coef.sig <- function(coef.names){
  # INPUTS:
  # coef.names: character string of coefficient names
  #
  # OUTPUTS:
  # name or index of last coefficient in table from lme output
  last.coef.name <- rev(coef.names)[1] # get name of last coefficient
  if(grepl('*',last.coef.name)){ # if last coefficient is an interaction
    last.coef.name <- gsub('\\*',':',last.coef.name)
  }
  return(last.coef.name)
}

lme.stepup <- function(mdl.gold=NULL,test.cfg=NULL,lme.ms.fun=lme.ms){
  # INPUTS:
  # mdl.gold: current gold standard model
  # test.cfg: list of fixed and random effects to step through
  # lme.ms: model function that allows formulas to be evaluated as plain text
  #
  # OUTPUTS:
  # gold standard model and associated config
  
  # start with first element of test config and only advance if that model is GOOD
  # application is to start with simple models and only step up if each model is better than the next
  test.counter <- 1 # start at first and advance as long as : (c1) p < 0.05 for anova(m1,m2) AND (c2) AIC drops AND (c3) added coefficient is significant
  cond <- TRUE
  while(cond & test.counter < length(test.cfg)+1){ # while cond is still true keep adding stuff until you've tested all models
    print('testing:')
    cfg <- test.cfg[[test.counter]]
    print(rev(cfg$fixed)[1])
    fitSuccess <- tryCatch({mdl.test <- lme.ms.fun(cfg,df = mdl.gold$data)
        cond <- lme.compare(mdl.gold,mdl.test,coef.sig=lme.last.coef.sig(cfg$fixed))
            },error = function(err) {
              print(err)
              print('halting model advancement now')
              return(FALSE)
            })
    if(fitSuccess){
      if(cond){
        mdl.gold <- mdl.test
        cfg.gold <- cfg} # if test successful update gold standard}
      test.counter <- test.counter+1
    } else {cond <- FALSE}
  }
  return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))
}

lme.stepdown <- function(mdl.gold=NULL,test.cfg=NULL,lme.ms.fun=lme.ms){
  # INPUTS:
  # mdl.gold: current gold standard model
  # test.cfg: list of fixed and random effects to step through
  # lme.ms: model function that allows formulas to be evaluated as plain text
  #
  # OUTPUTS:
  # gold standard model and associated config
  
  # start with first element of test config and only advance if that model is BAD
  # application is to start with complex models then step down until you find one better than gold
  test.counter <- 1 # start at first and advance as long as : NOT [ (c1) p < 0.05 for anova(m1,m2) AND (c2) AIC drops AND (c3) added coefficient is significant]
  adv <- TRUE # this tells you whether to keep advancing
  while(adv & test.counter < length(test.cfg)+1){ # while cond is still true keep adding stuff until you've tested all models
    print('testing:')
    cfg <- test.cfg[[test.counter]]
    print(rev(cfg$fixed)[1])
    fitSuccess <- tryCatch({mdl.test <- lme.ms.fun(cfg,df = mdl.gold$data)
    cond <- lme.compare(mdl.gold,mdl.test,coef.sig=lme.last.coef.sig(cfg$fixed)) # check if tested model is better than gold std
    },error = function(err) {
      print(err)
      print('halting model advancement now')
      return(FALSE)
    })
    if(fitSuccess){
      if(cond){ # if tested model IS better than gold, then store results and stop advancing
        mdl.gold <- mdl.test
        cfg.gold <- cfg # if test successful update gold standard
        adv <- FALSE
      } else if(!cond){adv <- TRUE} # if tested model is NOT better, then keep going
    } else {adv <- TRUE} # if you couldn't fit the first model, try again with the next less complex model
    test.counter <- test.counter+1
  }
  return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))
}

lme.selectbest <- function(mdl.gold=NULL,cfg.gold=NULL,test.cfg=NULL,lme.ms.fun=lme.ms){
  # INPUTS:
  # mdl.gold: current gold standard model
  # cfg.gold: configuration of current gold standard model
  # test.cfg: list of fixed and random effects to step through
  # lme.ms: model function that allows formulas to be evaluated as plain text
  #
  # OUTPUTS:
  # gold standard model and associated config

  mdl.tests <- lapply(test.cfg, function(cfg) tryCatch({return(lme.ms(cfg,df = mdl.gold$data))},
                                          error=function(err){
                                            return('')
                                          }))
  if(sum(mdl.tests != '')>0){ # if any of the tested models fit successfully
    mdl.tests <- mdl.tests[mdl.tests != ''] # delete models that couldn't be fit
    # first check if any models are doing better than gold standard
    test.vs.gold <- sapply(1:length(mdl.tests), function(mdl.idx) lme.compare(mdl.gold,mdl.tests[[mdl.idx]],coef.sig=lme.last.coef.sig(test.cfg[[mdl.idx]]$fixed)))
    #mdl.all <- c(list(mdl.gold),mdl.tests) # these should be in order
    # first column is whether each model is better than mdl.gold
    if(!any(test.vs.gold)){ # if no test model is better than mdl.gold, return mdl.gold and its cfg list
      return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))
    } else{ # otherwise, of the models better than gold, return the model that was better than the most other models tested 
      mdl.tests <- mdl.tests[test.vs.gold] # delete models not better than gold standard
      test.cfg <- test.cfg[test.vs.gold] # delete the configuration list elements associated with those models
      perf.comp <- matrix(NA,nrow=length(mdl.tests),ncol=length(mdl.tests)) # pairwise comparison of models that beat the gold standard
      for(mdl.idx1 in 1:length(mdl.tests)){
        for(mdl.idx2 in 1:length(mdl.tests)){
          perf.comp[mdl.idx1,mdl.idx2] <- lme.compare(mdl.tests[[mdl.idx2]],mdl.tests[[mdl.idx1]],coef.sig=lme.last.coef.sig(test.cfg[[mdl.idx1]]$fixed))
        }
      } # output mdl1 x mdl2 matrix of comparisons where rows indicate model being tested for coefficient significance and AIC drop and anova      
      mdl.idx <- which.max(rowSums(perf.comp))
      if(length(mdl.idx) > 1){print('ERROR: models tied'); return()}
      if(length(mdl.idx) == 0){print('ERROR: all tested models are equivalent'); return()}
      mdl.gold <- mdl.tests[[mdl.idx]]
      cfg.gold <- test.cfg[[mdl.idx]]
      return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))
    }
  } else {return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))}
}


lme.fixed.predict.ave <- function(mdl,fixed.terms){
  # INPUTS:
  # mdl: lme model object
  # fixed.terms: characters specifying which terms to get fixed effect overall line (no random effects)
  # using group average values for all other covariates
  #
  # OUTPUTS:
  # for every combination of levels in fixed.terms
  # get fitted values of model using fixed effects and average values for covariates 
  # return response variable averaged as well as the fitted values with .fixed appended
  
  df <- mdl$data[complete.cases(mdl$data),] # get complete cases
  df.agg <- aggregate(df,by=setNames(lapply(fixed.terms, function(term) df[,term]),fixed.terms),mean)
  df.agg[,'(Intercept)'] <- 1 # add intercept
  fixed.coefs <- mdl$coefficients$fixed
  fixed.coef.names <- names(fixed.coefs)
  if(any(grepl(':',fixed.coef.names))){ # if there are any interactions, make a new column in df.agg for the interaction
    for(int.idx in grep(':',fixed.coef.names)){
      coef.int <- fixed.coef.names[int.idx] # get name of interaction (var1:var2)
      int.components <- strsplit(coef.int,':',fixed = T)[[1]] # split at colon
      # create interaction predictor column by multiplying variables together
      # doing this with averaging only really makes sense with reasonably discrete variables
      df.agg[,coef.int] <- df.agg[,int.components[1]]*df.agg[,int.components[2]]
    }
  }
  resp.name <- attr(getResponse(mdl),'label')
  df.agg[,paste0(resp.name,'.fixed')] <- sapply(1:nrow(df.agg), function(j) sum(fixed.coefs * df.agg[j,fixed.coef.names]))
  return(df.agg[,c(resp.name,paste0(resp.name,'.fixed'),fixed.terms)])
}

lme.subj.coefs <- function(mdl,c.o.i){
  # INPUTS:
  # mdl: lme model object
  # c.o.i: character specifying which terms to get subject-specific parameters for
  #
  # OUTPUTS:
  # for coefficients in c.o.i
  # get coefficients for each subject by combining fixed and random effects, and interactions if present
  # return dataframe of subject-specific coefficients with subject identifiers
  # ** assumes you only have one id variable for random effects **
  
  fixed.coefs <- mdl$coefficients$fixed
  random.coefs <- mdl$coefficients$random[[1]] # ** assumes you only have one id variable for random effects **
  subj.coefs <- intersect(names(fixed.coefs),colnames(random.coefs)) 
  return(fixed.coefs[subj.coefs] + random.coefs[2,subj.coefs])
  
  if(any(grepl(':',fixed.coef.names))){ # if there are any interactions, make a new column in df.agg for the interaction
    for(int.idx in grep(':',fixed.coef.names)){
      coef.int <- fixed.coef.names[int.idx] # get name of interaction (var1:var2)
      int.components <- strsplit(coef.int,':',fixed = T)[[1]] # split at colon
      # create interaction predictor column by multiplying variables together
      # doing this with averaging only really makes sense with reasonably discrete variables
      df.agg[,coef.int] <- df.agg[,int.components[1]]*df.agg[,int.components[2]]
    }
  }
  resp.name <- attr(getResponse(mdl),'label')
  df.agg[,paste0(resp.name,'.fixed')] <- sapply(1:nrow(df.agg), function(j) sum(fixed.coefs * df.agg[j,fixed.coef.names]))
  return(df.agg[,c(resp.name,paste0(resp.name,'.fixed'),fixed.terms)])
  
  
}
