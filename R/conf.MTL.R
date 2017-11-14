#' Confidence interval of mean total lifetime
#'
#' This function is used to calculate confidence intervals of mean total lifetime using jackknife resampling.
#'
#' This function systematically leaves out each subject from the original datset and simulates mean total lifetimes
#' for each \code{n-1}-sized subsample. The jackknife mean and variance are calculated by aggregating \code{n} simulated
#' mean total lifetimes. For each jackknife dataset, mean total lifetime is simulated using the
#' algorithm described in sim.MTL.
#'
#'
#' @param obj An object returned by \code{optim.fit}, which contains the transition probabilities
#' and other information used to simulate mean total lifetime.
#' @param state A numeric vector indicating from which state the mean total lifetime is simulated.
#' Default is NULL, where no mean total life for a specific state is output. If obj is returned by optim.fit with
#' treatment=NULL, there is no need to set this argument.
#' @param nsim The times of simulation for mean total life. The default is 1000.
#' @param L The prespecified threshold for blocking the increase of residual lifetime. The default is 120.
#' @return If the input object comes from \code{optim.fit} with \code{treatment=NULL}, a list object with elements:
#' \item{conf.state.MTL}{A data frame containing states, corresponding mean total lifetime,  standard
#'  error and 95\% confidence interval. If state=NULL, this element does not exist.}
#' \item{state.table}{The correspondence of state number and state label.}
#' If the input object comes from \code{optim.fit} with \code{treatment} is not NULL, a list object with elements:
#' \item{conf.strategies}{Mean total lifetime for different strategies, along with standard error and 95\% confidence interval}
#'
#' @aliases conf.MTL
#' @seealso \code{\link{optim.fit}}
#' @export
#' @examples
#'
#'\dontrun{
#' library(OptimalTiming)
#'
#'##################################
#'## Example 1: This example shows how to calculate confidence
#'## intervals for different treatment strategies
#'
#'## read data
#'data(SimCml)
#'
#'## fit multistate model with treatment not equals NULL
#'fit=optim.fit(data=SimCml,
#'        transM=matrix(c(0,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,1,1,1,
#'        0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0),7,byrow=TRUE),
#'        nstate=7,state_label=c("diagnose","cp1","ap","cp2","bc","sct","death"),
#'        event_label=c("cp1.s","ap.s","cp2.s","bc.s","sct.s","death.s"),
#'        treatment=c("sct","sct.s"),absorb=c("death","death.s"),
#'        cov=c("age"),cov_value=c(0))
#'
#'## compare different treatment strategies
#'conf.MTL(obj=fit,nsim=1000,L=120)
#'
#'##################################
#'## Example 2: This example shows how to calculate confidence
#'## intervals for a given state
#'
#'## read data
#'data(SimCml)
#'
#'## delete the information of transplant time
#'data=SimCml[SimCml$sct.s==0,]
#'del=which(names(SimCml)%in%c("sct","sct.s"))
#'data=data[,-del]
#'
#'## fit multistate model with treatment equals NULL
#'fit=optim.fit(data=data,
#'         transM=matrix(c(0,1,0,0,0,0,0,0,1,0,1,1,0,0,0,
#'         1,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0),6,byrow=TRUE),
#'         nstate=6,state_label=c("diagnose","cp1","ap","cp2","bc","death"),
#'         absorb=c("death","death.s"),event_label=c("cp1.s","ap.s","cp2.s","bc.s","death.s"),
#'         cov=c("age"),cov_value=c(0))
#'
#'## calculate mean total lifetime and confidence intervals
#'## for state 1,2,3,4
#'conf.MTL(obj=fit,state=c(1,2,3,4),nsim=1000,L=120)}
#'
#'
#'



conf.MTL <- function(obj,state=NULL,nsim=1000,L=120){
##########################################################
## treatment = NULL
  if(is.null(obj$overall$treatment))
  {
    data=obj$overall$data
    transM=obj$overall$transM
    nstate=obj$overall$nstate
    state_label=obj$overall$state_label
    event_label=obj$overall$event_label
    cov=obj$overall$cov
    cov_value=obj$overall$cov_value

    options(warn=1)
    jacknife=dim(data)[1]
    conf_result.state=vector("list",jacknife)
    if(is.null(state))
      {
      conf_result.state[[i]]=NULL
      }else
      {
        cat("Jackknife sampling will take several hours, please wait...\n")
        for(i in 1:jacknife){
        bs_data=data[-i,]
        ofit=optim.fit(data=bs_data,transM=transM,nstate=nstate,state_label=state_label,event_label=event_label,treatment=NULL,absorb=obj$overall$absorb,cov=cov,cov_value=cov_value)
        tobj=sim.MTL(ofit,state=state,nsim=nsim,L=L)
        conf_result.state[[i]]=tobj$state.MTL
        print(i)
        }
      }

    ## conf.state.MTL
    if(!is.null(state))
    {
      bs_result=do.call(rbind,conf_result.state)
      bs_result=bs_result[order(bs_result[,1]),]
      bs_mean_state=as.matrix(tapply(bs_result[,2],bs_result[,1],mean))
      bs_sd_state=as.matrix(tapply(bs_result[,2],bs_result[,1],sd))
      upci=bs_mean_state+bs_sd_state*qnorm(0.975)
      lowci=bs_mean_state-bs_sd_state*qnorm(0.975)
      conf.state.MTL=data.frame(mean=bs_mean_state,sd=bs_sd_state,lowCI=lowci,upCI=upci)
    }

    ## output
    if(is.null(state)){
      object=list()
      object
    }else
    {
      object=list()
      object$conf.state.MTL=conf.state.MTL
      object$state.table=data.frame(state=1:obj$overall$nstate,label=obj$overall$state_label)
      object
    }

##########################################################
## treatment != NULL
  }else if (!is.null(obj$overall$treatment))
  {
    data=obj$overall$data
    transM=obj$overall$transM
    nstate=obj$overall$nstate
    state_label=obj$overall$state_label
    event_label=obj$overall$event_label
    cov=obj$overall$cov
    cov_value=obj$overall$cov_value

    options(warn=1)
    jacknife=dim(data)[1]
    conf_result.strategies=vector("list",jacknife)

    cat("Jackknife sampling will take several hours, please wait...\n")
    for(i in 1:jacknife)
    {
      bs_data=data[-i,]
      ofit=optim.fit(data=bs_data,transM=transM,nstate=nstate,state_label=state_label,event_label=event_label,treatment=obj$overall$treatment,absorb=obj$overall$absorb,cov=cov,cov_value=cov_value)
      conf.obj=sim.MTL(ofit,nsim=nsim,L=L)
      conf_result.strategies[[i]]=conf.obj$strategies
      print(i)

    }

    ## conf.strategies
    bs_result=do.call(rbind,conf_result.strategies)
    bs_mean=as.matrix(tapply(bs_result[,2],bs_result[,1],mean))
    bs_sd=as.matrix(tapply(bs_result[,2],bs_result[,1],sd))
    upci=bs_mean+bs_sd*qnorm(0.975)
    lowci=bs_mean-bs_sd*qnorm(0.975)
    conf.strategies=data.frame(mean=bs_mean,sd=bs_sd,lowCI=lowci,upCI=upci)

    object=list()
    object$conf.strategies=conf.strategies
    object
  }



}

