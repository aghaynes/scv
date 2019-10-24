#' Benchmark funnel plot boundaries
#'
#' This function produces data for drawing funnel plot boundaries that correspond to the 95\% confidence intervals
#' of a Poisson-gamma distribution corresponding to a desired systematic component of variation, the benchmark SCV.
#'
#' @param exp expected values from indirect standardization (observed/expected)
#' @param scvtx benchmark SCV , by default 3
#' @param alphax alpha level, by default 95\%
#'
#' @return exp
#'
#' upperCI Upper limit O/E of benchmark Poisson-Gamma distribution
#'
#' lowerCI Lower limit O/E of benchmark Poisson-Gamma distribution
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
#' outlier<-outlier_stats(obs,exp,scvtx=3, alphax=.95)
#'
#' funnel<-benchmark_funnel(exp,scvtx=3, alphax=.95)
#'
#' ggplot(funnel,aes(exp,upperCI))+
#' geom_line()+
#' geom_line(aes(y=lowerCI))+
#' geom_point(data=outlier$outlier_stats,aes(x=exp, y=obs/exp,color=is_outlier))+
#' xlab("Expected cases")+
#' ylab("Observed/Expected")+
#' labs(color="Outlier")+
#' theme_minimal()
#'
#' @export
#' @import DCluster ggplot2

benchmark_funnel<- function(exp,scvtx=3,alphax=.95){

  step<-ifelse(max(exp)<50,1,
               ifelse(max(exp)<200,5,10))

  px<-(1-alphax)/2

  nbinter<-(100/scvtx)/ ((100/scvtx)+seq(1,max(exp)+1,step))


  qu1<-qnbinom(p=1-px,size = (100/scvtx),prob = nbinter)
  pqu1<-pnbinom(qu1,size =(100/scvtx),prob = nbinter)
  qu2<-qnbinom(p=1-px,size = (100/scvtx),prob = nbinter)-1
  pqu2<-pnbinom(qu2,size = (100/scvtx),prob = nbinter)
  a<-((1-px) - pqu2)/(pqu1-pqu2)
  uppernegbin<-(qu2 + a*(qu1-qu2))/seq(1,max(exp)+1,step)



  ql1 <- qnbinom(px, size = (100/scvtx),prob = nbinter) - 1
  l1 <- replace(ql1, ql1<0, 0)
  p.l1 <- pnbinom(l1, size = (100/scvtx),prob = nbinter)
  l2 <- qnbinom(px,  size = (100/scvtx),prob = nbinter)
  p.l2 <- pnbinom(l2, size = (100/scvtx),prob = nbinter)
  a <- (px - p.l1)/(p.l2-p.l1)
  lowernegbin <- (l1 + a*(l2-l1))/seq(1,max(exp)+1,step)
  lowernegbin<-replace(lowernegbin,is.na(lowernegbin),0)

  ic<-cbind("exp"=seq(1,max(exp)+1,step),uppernegbin,lowernegbin)
  colnames(ic)<-c("exp","upperCI","lowerCI")

  as.data.frame(ic)}


