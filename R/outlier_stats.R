#' Outlier statistics and outlier gains
#'
#' Outlier detection from  benchmark Poisson-Gamma distribution (a distribution drawn from a benchmark SCV).
#' Outlier gains are the additions or reductions (depending on whether it is a desired outcome or not) to the system that would occur
#' under a scenario where all identified upper outliers or lower outliers are brought into the 95\% bounds of the benchmark distribution
#'
#'
#' @param obsx  observed cases for a given population
#' @param expx expected cases from standardization of given population
#' @param scvtx benchmark SCV, default 3
#' @param alphax alpha level, default .95
#'
#'
#' @import DCluster ggplot2
#'
#' @return obs
#'
#' exp
#'
#' upper_lim_qnbinom Upper limit O/E of benchmark Poisson-Gamma distribution
#'
#' lower_lim_qnbinom Lower limit O/E of benchmark Poisson-Gamma distribution
#'
#' is_outlier Identifes O/E as upper or lower outlier (else is NULL)
#'
#' pval_ref Probability of O/E ratio in benchmark Poisson-Gamma distribution
#'
#' outlier_diff Difference between outlier and benchmark upper or lower limit
#'
#' upper_outlier_gain System gains under scenario of normalizing upper outliers
#'
#' lower_outlier_gain System gains under scenario of normalizing lower outliers
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
#' outlier_stats(obs,exp,scvtx=3, alphax=.95)
#'
#' @export


outlier_stats <-function(obsx,expx, scvtx=3, alphax=.95){

  list(outlier_stats=data.frame(
    "obs"=obsx,
    "exp"=expx,
    "upper_lim_qnbinom"=.qnbinom_lim_func(expx,scvtx,alphax)$qnbinom_upp/expx,
    "lower_lim_qnbinom"=.qnbinom_lim_func(expx,scvtx,alphax)$qnbinom_low/expx,
    "is_outlier" = as.character(.outlier_func_ref(obsx,expx,scvtx,alphax)),
    "pval_ref"=.pvf_ref(obsx,expx,scvtx),
    "outlier_diff"=.outlier_diff_func(obsx,expx,scvtx,alphax),

    stringsAsFactors = FALSE),
    outlier_gains=c( "upper_outlier_gain" = .upper_outlier_gain_func(obsx,expx,scvtx,alphax),
                     "lower_outlier_gain" = .lower_outlier_gain_func(obsx,expx,scvtx,alphax)))

}

