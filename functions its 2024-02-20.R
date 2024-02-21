#Adapted from Weinberger Vaccine eval workskshop, Gambia example
#Instead of useing three post-vaccine periods, only using one post-COVID period

step_func <- function(ds,
                     time_points=unique(ds$date),
                     post_period1=c('2009-04-01', '2011-12-31'),
                     post_period2=c('2012-01-01', '2019-12-31'),
                     vax.vars=c('post1', 'post2', 'post1.time', 'post2.time'),
                     other.covars,
                     mod,
                     denom, 
                     outcome_name){
  
  ds <- ds[order(ds$date),]
  ds$time_index<-(1:nrow(ds))/nrow(ds)
  
  #Create the dummy variable for PCV period
  ds$post1<-as.numeric(ds$date >= post_period1[1] & ds$date <  post_period1[2]) 
  ds$post2<-as.numeric(ds$date >= post_period2[1] & ds$date <  post_period2[2]) 
  
  ds$time_index<-1:nrow(ds)
  ds$sin12<-sin(2*pi* ds$time_index/12)
  ds$cos12<-cos(2*pi* ds$time_index/12)
  ds$sin6<-sin(2*pi* ds$time_index/6)
  ds$cos6<-cos(2*pi* ds$time_index/6)
  ds$sin3<-sin(2*pi* ds$time_index/3)
  ds$cos3<-cos(2*pi* ds$time_index/3)
  ds$month<-month(ds$date)
  ds$time_index<-ds$time_index/max(ds$time_index)
  
  #Time index variables
  ds$post1.time <- ds$post1*(ds$time_index)
  ds$post2.time <- ds$post2*(ds$time_index)
  
    #Create month dummies if scale<annual
  if((ds %>% 
      mutate(year = year(date)) %>% 
      group_by(year) %>% 
      summarise(n=n()))$n[1]>1){
    month.dummies <- model.matrix(~ as.factor(month), data=ds)
      month.dummies <- as.data.frame(month.dummies[,-1])
      names(month.dummies) <- paste0('month', 2:(ncol(month.dummies)+1))
      ds<-cbind.data.frame(ds, month.dummies)
      
      seas.vars<- c(names(month.dummies) )
  }else{seas.vars<- c()}
 
  ds$one<-1
  ds$obs <- as.factor(1:nrow(ds))
  #ds$log.offset<-log(ds[,denom]+0.5)
  #ds$log.offset<-log(ds$pop+0.5)
  ds <- ds %>% 
    mutate(log.offset = -log(!!sym(denom)+0.5))
 #ds$log.offset<-log(ds$!!sym(denom)+0.5)
  
  offset.vars<-'offset(log.offset)'
  
  if((other.covars=='none')[1]){
    mod.vars<-c(seas.vars,vax.vars)
    mod.covars.cf <- seas.vars
  }else{
    mod.vars<-c(seas.vars,other.covars,vax.vars)
    mod.covars.cf<-c(seas.vars,other.covars)
    #for(i in 1:length(other.covars)){
    #  ds[,other.covars[i]]<-scale(log( ds[,other.covars[i]]+0.5)) 
    #}
  }
  
  #Re-define time_index (it gets overwritten in loop above if other covars present)
  ds$time_index<-(1:nrow(ds))/nrow(ds)
  
  form1<-as.formula(paste0(outcome_name,'~',paste0(c(mod.vars,offset.vars), collapse='+')))
  
  #Fit model
  if(mod=='pois'){
  mod1 <-
    glm(form1,
        data = ds, family='poisson'
        )
  }else if(mod=='negbin'){
    mod1 <-
      glm.nb(form1,
          data = ds
      )
  }

  aic1=mod1$aic[[1]]
  deviance1=mod1$deviance[[1]]
  df1=mod1$df.residual[[1]]
  overdispersion=deviance1/df1

  #GENERATE PREDICTIONS
  covars3 <-
    as.matrix(cbind(ds[, mod.vars])) 
  covars3 <- cbind.data.frame(rep(1, times = nrow(covars3)), covars3)
  names(covars3)[1] <- "Intercept"
  
  pred.coefs.reg.mean <-
    mvrnorm(n = 1000,
            mu = coef(mod1),
            Sigma = vcov(mod1))
  
  preds.stage1.regmean <-
    exp(as.matrix(covars3) %*% t(pred.coefs.reg.mean) +ds$log.offset)
  
  preds.q<-t(apply(preds.stage1.regmean,1,quantile, probs=c(0.025,0.5,0.975)))
  
  #Then for counterfactual, set post-vax effects to 0.
  covars3.cf <-
    as.matrix(cbind(ds[, c(
      mod.covars.cf
    )], matrix(
      0, nrow = nrow(ds), ncol = length(vax.vars)
    )))
  covars3.cf <-
    cbind.data.frame(rep(1, times = nrow(covars3.cf)), covars3.cf)
  preds.stage1.regmean.cf <-    exp(as.matrix(covars3.cf) %*% t(pred.coefs.reg.mean)+ ds$log.offset)
  preds.cf.q<-t(apply(preds.stage1.regmean.cf,1,quantile, probs=c(0.025,0.5,0.975)))
  
  rr.t <- preds.stage1.regmean / preds.stage1.regmean.cf
  rr.q.t <- t(apply(rr.t, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  
  last.t<-nrow(rr.t) #evaluate at last time point
  preds.stage1.regmean.SUM <-   preds.stage1.regmean[last.t, ]
  preds.stage1.regmean.cf.SUM <-preds.stage1.regmean.cf[last.t, ]
  rr.post <- preds.stage1.regmean.SUM / preds.stage1.regmean.cf.SUM
  rr.q.post <- quantile(rr.post, probs = c(0.025, 0.5, 0.975))
  
  #Combined dataset to plot from
  mod1.ds <- cbind(ds, preds.q[,'50%'], preds.cf.q[,'50%'])
  
  #For yearly RRs for Years
  preds.cf.annual <-  cbind.data.frame(date=year(unique(ds$date)),preds.cf.q) %>% 
    group_by(date) %>% 
    summarise(preds.cf.q.2.5 = sum(`2.5%`),
              preds.cf.q.50 = sum(`50%`),
              preds.cf.q.97.5 = sum(`97.5%`))
  
  preds.annual <-  cbind.data.frame(date=year(unique(ds$date)),preds.q) %>% 
    group_by(date) %>% 
    summarise(preds.q.2.5 = sum(`2.5%`),
              preds.q.50 = sum(`50%`),
              preds.q.97.5 = sum(`97.5%`))
  
  #Denominators
  denom <- ds %>% 
    group_by(year(date)) %>% 
    summarise(denom = max(pop))
  
  ipd_true <- ds %>% 
    group_by(year(date)) %>% 
    summarise(ipd_true = sum(ipd))
  
 rr.annual <- preds.annual/preds.cf.annual
 names(rr.annual) <- c("obsv", "rr.q.2.5", "rr.q.50", "rr.q.97.5")
 rr.annual$date <- as.numeric(rownames(rr.annual))
 rr.annual$ds <- as.character(k)
 rr.annual <- bind_cols(rr.annual, preds.annual, preds.cf.annual, denom, ipd_true)
 
 #Add rate
 #rr.annual <- rr.annual %>% 
#   mutate(ipd_inc_true = (ipd_true/denom)*100000,
#          ipd_inc_exp = (rr.annual$preds.cf.q.50/rr.annual$denom)*100000)
   
  
 rr.out <- list('rr.q.post' = rr.q.post, 'outcome'=ds[,outcome_name],
                'preds.cf.q'=preds.cf.q,'preds.q'=preds.q,
                'rr.q.t'=rr.q.t,'overdispersion'=overdispersion, 'dates'=ds$date, 'mod1.ds'=mod1.ds, "rr.annual"=rr.annual,
                'aic' = aic1, 'deviance' = deviance1, 'df' = df1, 'overdispersion' = overdispersion)
 rr.out
  
}

#Functions for plotting
my_ceil <- function(x) {
  ceil <- ceiling(max(x))
  ifelse(ceil > 1 & ceil %% 2 == 1, ceil + 1, ceil)
}

my_breaks <- function(x) {
  ceil <- my_ceil(max(x))
  if (ceil < 2) {
    breaks <- seq(0, ceil, by = 0.5)
  } else {
    breaks <- ceiling(pretty(seq(0, ceil)))
  }
  unique(breaks)
}

my_limits <- function(x) { 
  ceil <- my_ceil(x[2])
  c(x[1], ceil)
}
