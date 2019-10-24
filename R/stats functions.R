#' @title Internal functions
#' @description Internal functions
#'
#' @keywords internal
#' @param obsx observed
#' @param expx expected
#' @param scvtx benchmark scv
#' @param alphax threshold level
#'
#' @keywords internal
#'
#' @import DCluster ggplot2
#'



### #systematic variance from poisson gamma model
.pg_var<-function(obsx,expx){
  mp<-empbaysmooth(obsx,expx,maxiter=200)
  mp$nu/mp$alpha^2
 }


#simulation set
.simf<-function(obsx, expx){
  nsim<-1000
  lambda0 <- sum(obsx)/sum(expx)

  ebs<-empbaysmooth(obsx,expx,maxiter=200)
  alpha<-ebs$alpha
  nu<-ebs$nu


  sim.lambda <- matrix(rgamma(length(expx) * nsim, shape = nu, scale = 1/alpha ),
                       nrow = length(expx))

  sim.z <- matrix(rpois(nrow(sim.lambda) * nsim, expx *
                          sim.lambda), nrow = length(expx))


  ebsim<-sapply(1:nsim,function(i){
    if(sum(sim.z[, i])>0 &&
       var(sim.z[, i][expx>0]/expx[expx>0])!=0){
      if(!is.na(tryCatch(empbaysmooth(sim.z[,i],expx,maxiter=200)[[1]],
                         error=function(err) NA))){
        out <-empbaysmooth(sim.z[, i], expx,maxiter=200)
        data.frame("var"=out$nu/out$alpha^2)

        }else{data.frame("var"=NA)
      }}else{data.frame("var"=NA)
    }
  },simplify =TRUE)
 ebsim
}

#variance of the scv
.var_var<-function(simx){
  us<-unlist(simx)
  if(length(us[is.na(us)])/length(us)<.01){
  var(us[!is.na(us)])}else{NA}
}



#Bootstrapp scv 95% ci
.pg_var_ci<-function(simx){
  us<-unlist(simx)
  if(length(us[is.na(us)])/length(us)<.01){
    list("lower"=100*as.numeric(quantile(us[!is.na(us)],probs=.025)),
         "upper"=100*as.numeric(quantile(us[!is.na(us)],probs=.975)))
  }else{list("lower"=NA,
             "upper"=NA)}
  }



# #probability from empirical distribution
.pvf<-function(obsx,expx){
  es<-empbaysmooth(obsx,expx,maxiter = 200)
  probx<-es$alpha/ (es$alpha+expx)
  pnbinom(q=obsx,size =es$nu,prob = probx)
}

#probability from benchmark distribution
.pvf_ref<-function(obsx,expx,scvtx){
  probx<-(100/scvtx)/ ((100/scvtx)+expx)
  pnbinom(q=obsx,size =(100/scvtx),prob = probx)
}


# #outlier from benchmark dist
.outlier_func_ref<-function(obsx,expx,scvtx,alphax){
  px<-(1-alphax)/2
  ifelse(obsx/expx > 1,
         ifelse(.pvf_ref(obsx,expx,scvtx = scvtx)>=(1-px),"upper","NULL"),
         ifelse(.pvf_ref(obsx,expx,scvtx = scvtx)<=px,"lower","NULL"))
}


# #smooth ratios
.smooth_func<- function(obsx,expx){
  empbaysmooth(obsx,expx,maxiter = 200)$smthrr}




# #smoothed observed
.obs_smooth_func<-function(obsx,expx){expx*.smooth_func(obsx,expx)}


#limits of qnbinom from reference distribution
.qnbinom_lim_func<-function(expx,scvtx,alphax){
  px<-(1-alphax)/2
  probx<-(100/scvtx)/ ((100/scvtx)+expx)
  qnbinom_upp<-qnbinom(p=1-px,size =(100/scvtx),prob = probx)
  qnbinom_low<-qnbinom(p=px,size =(100/scvtx),prob = probx)
  list("qnbinom_upp"=qnbinom_upp,"qnbinom_low"=qnbinom_low)
}


#outlier difference based on observed and qnbinom
.outlier_diff_func<-function(obsx,expx,scvtx,alphax){
  px<-(1-alphax)/2
  probx<-(100/scvtx)/ ((100/scvtx)+expx)
  qnbinom_upp<-qnbinom(p=1-px,size =(100/scvtx),prob = probx)
  qnbinom_low<-qnbinom(p=px,size =(100/scvtx),prob = probx)
  ifelse(qnbinom_upp-obsx<0, qnbinom_upp-obsx,
         ifelse(qnbinom_low-obsx>0,qnbinom_low-obsx,0))
}


# outlier gains (gains from bringing outliers to upper Confidence interval boundary of reference distribution)
.upper_outlier_gain_func<-function(obsx,expx,scvtx,alphax){
 outlierdiffs<-.outlier_diff_func(obsx,expx,scvtx,alphax)
 sum(outlierdiffs[outlierdiffs<0])
}


# inferieur outlier gains (gains from bringing smoothed observed to lower Confidence interval boundary of reference distribution)
.lower_outlier_gain_func<-function(obsx,expx,scvtx,alphax){
  outlierdiffs<-.outlier_diff_func(obsx,expx,scvtx,alphax)
  sum(outlierdiffs[outlierdiffs>0])
}

#







