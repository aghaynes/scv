#' Poisson-Gamma systematic component of variation 
#'
#'  
#' @param obsx  observed
#' @param expx expected
#' 
#' @import DCluster
#' 
#' @return obs
#'
#'  exp
#' 
#' pval probability of O/E ratio in empirical Poisson-Gamma distribution
#' 
#' smoothed empirical bayes smoothing of O/E ratio
#' 
#' scv systematic component of variation
#' 
#' scv_var bootstrapped variance of the scv
#' 
#' scv_CI_low bootstrapped lower 95% confidence interval of SCV
#' 
#' scv_CI_upp bootstrapped upper 95% confidence interval of SCV
#' 
#' 
#' @examples
#'
#' data(nc.sids)
#' 
#' obs<-nc.sids$SID74
#' 
#' exp<-nc.sids$BIR74*sum(nc.sids$SID74)/sum(nc.sids$BIR74)
#' 
#' scv_stats(obs,exp)
#' 
#' 
#' @export


scv_stats<-function(obsx,expx){
 
  sim_set<-suppressWarnings(.simf(obsx,expx))
  
  list("smoothing"=data.frame(
    "obs"=obsx,
    "exp"=expx,
    "pval" =  .pvf(obsx,expx),
    "smoothed"=.smooth_func(obsx,expx)),
    "scv_stats"=
      c("scv" = 100*.pg_var(obsx,expx),
        "scv_var"=.var_var(sim_set),
        "scv_CI_low"=.pg_var_ci(sim_set)$lower,
        "scv_CI_upp"=.pg_var_ci(sim_set)$upper))
}

