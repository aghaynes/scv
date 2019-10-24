#' Outlier detection from empirical and benchmark Poisson-Gamma distribution
#'
#'
#' @param obsx  observed
#' @param expx expected
#' @param scvtx benchmark SCV, default 3
#' @param alphax alpha level, default .95
#'
#'
#' @import DCluster
#'
#' @return obs
#'
#' exp
#'
#' upper_lim_qnbinom upper limit O/E of benchmark Poisson-Gamma distribution
#'
#' lower_lim_qnbinom lower limit O/E of benchmark Poisson-Gamma distribution
#'
#' is_outlier identifes O/E as upper or lower outlier (else is NULL)
#'
#' pval_ref probability of O/E ratio in benchmark Poisson-Gamma distribution
#'
#' outlier_diff difference between outlier and benchmark upper or lower limit
#'
#' upper_outlier_gain system gains under from scenario of normalizing upper outliers
#'
#' lower_outlier_gain system gains under from scenario of normalizing lower outliers
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

