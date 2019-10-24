#' Poisson-Gamma systematic component of variation
#'
#' Extension of the empbaysmooth function in DCluster. Adds a computation of the SCV as 100* \emph{nu}/\emph{alpha^2}.
#' Computes from simulations a bootstrapped variance and 95\% CI of the SCV.
#'
#' @param obsx  observed cases for a given population
#'
#' @param expx expected cases from standardization of given population
#'
#' @import DCluster ggplot2
#'
#' @return obs
#'
#'  exp
#'
#' pval Probability of O/E ratio in empirical Poisson-Gamma distribution
#'
#' smoothed Empirical bayes smoothing of O/E ratio
#'
#' scv Systematic component of variation
#'
#' scv_var Bootstrapped variance of the scv
#'
#' scv_CI_low Bootstrapped lower 95\% confidence interval of SCV
#'
#' scv_CI_upp Bootstrapped upper 95\% confidence interval of SCV
#'
#'
#' @examples
#'
#' data(nc.sids, package= "spData")
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

