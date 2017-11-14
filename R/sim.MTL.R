#' Simulate mean total lifetime
#'
#' This function is used to simulate mean total lifetime for a given initial state
#' according to the estimated transition probabilities.
#'
#'
#' This part describes the algorithm used to simulate mean total lifetime in detail.
#' For an initial state, we first extract the transition probability data frame for this state,
#' and cumulate probabilities in row direction. In the transformed data frame, the first column
#' of which impels the time points where transition probabilities are measured, and the last column
#' is accumulated to be 1, indicating the addition of the chance to depart from a state and the
#' chance to remain at the state equals 1. Then, we generate a random
#' value from uniform distribution \code{Unif(0,1)} to determine the next state by comparing
#' this value to cumulative probabilities. The state, whose cumulative transition probability
#' from initial state firstly surpasses the uniform value, is defined as the next state.
#' The time interval, from the initiation of the study to where the transition take place,
#' is defined as the interim residual life. These two variables (next state and interim residual life)
#' are recorded for late use. Subsequently, we regard the next state as the
#' initial state, and repeat this searching process until the absorbing state is
#' reached or the interim residual lifetime surpasses the prespecified threshold (\code{L}).
#' Finally, the mean total life is either the last interim residual lifetime from the initiation
#' of study to the occurrence of absorbing state, or the prespecified threshold (\code{L}) if the
#' absorbing has not reached yet.
#'
#' If a state is given in this function, we set this state as initial state and perform the algorithm
#' mentioned above for \code{nsim} times, and average the output to obtain mean total lifetime.
#'
#'
#' According to different type of input object, this function return different results.
#' If the object comes from optim.fit with \code{treatment=NULL}, this function is used to simulate mean
#' total lifetime for a given state. If the object comes from optim.fit with \code{treatment} not
#' equals \code{NULL}, this function is used to compare mean total lifetimes of subjects who receive the new treatment to
#' those who do not receive the new treatment.
#'
#' @param obj An object returned by \code{optim.fit}, which contains the transition probabilities
#' and other information used to simulate mean total lifetime.
#' @param state A numeric vector indicating from which state the mean total lifetime is simulated.
#' Default is NULL, where no mean total life for a specific state is output. If obj is returned by optim.fit with
#' treatment=NULL, there is no need to set this argument.
#' @param nsim The times of simulation. The default is 1000.
#' @param L The prespecified threshold for blocking the increase of residual lifetime. The default is 120.
#' See Details.
#' @return
#' If the input object comes from \code{optim.fit} with \code{treatment=NULL}, a list object with elements:
#' \item{state.MTL}{A data frame containing states and corresponding mean total lifetime. If state=NULL, this element does not exist.}
#' \item{state.table}{The correspondence of state number and state label.}
#' If the input object comes from \code{optim.fit} with \code{treatment} not equals NULL, a list object with elements:
#' \item{strategies}{Mean total lifetime for different strategies.}
#'
#' @aliases sim.MTL
#' @seealso \code{\link{optim.fit}}
#' @export
#' @examples
#'
#'\dontrun{
#'library(OptimalTiming)
#'
#' ##################################
#' ## Example 1: This example shows how to use this package to find
#' ## the optimal timing of new treatment initiation
#'
#' ## read data
#' data(SimCml)
#'
#' ## fit multistate model with treatment not equals NULL
#' fit=optim.fit(data=SimCml,
#'          transM=matrix(c(0,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,1,1,1,0,0
#'          ,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0),7,byrow=TRUE),
#'          nstate=7,state_label=c("diagnose","cp1","ap","cp2","bc","sct","death"),
#'          event_label=c("cp1.s","ap.s","cp2.s","bc.s","sct.s","death.s"),
#'          treatment=c("sct","sct.s"),absorb=c("death","death.s"),
#'          cov=c("age"),cov_value=c(0))
#'
#' ## compare different treatment strategies
#' sim.MTL(obj=fit,nsim=1000,L=120)
#'
#'##################################
#' ## Example 2: This example shows how to obtain mean total lifetime
#' ## for a given state
#'
#' ## read data
#' data(SimCml)
#'
#' ## delete the information of transplant time
#' data=SimCml[SimCml$sct.s==0,]
#' del=which(names(SimCml)%in%c("sct","sct.s"))
#' data=data[,-del]
#'
#' ## fit multistate model with treatment equals NULL
#' fit=optim.fit(data=data,
#'         transM=matrix(c(0,1,0,0,0,0,0,0,1,0,1,1,0,0,0,1,
#'         1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0),6,byrow=TRUE),
#'         nstate=6,state_label=c("diagnose","cp1","ap","cp2","bc","death"),
#'         absorb=c("death","death.s"),
#'         event_label=c("cp1.s","ap.s","cp2.s","bc.s","death.s"),
#'         cov=c("age"),cov_value=c(0))
#'
#' ## calculate mean total lifetime when the initiate state is cp1 or ap
#' sim.MTL(obj=fit,state=c(2,3),nsim=1000,L=120)}
#'
#'\dontshow{
#'library(OptimalTiming)
#'data(SimCml)
#'SimCml=SimCml[SimCml$sct.s==0,]
#'SimCml=SimCml[,c(1,2,3,4,11,12,13,14)]
#'fit=optim.fit(data=SimCml,
#'          transM=matrix(c(0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0),4,byrow=TRUE),
#'          nstate=4,state_label=c("diagnose","cp1","ap","death"),
#'          event_label=c("cp1.s","ap.s","death.s"),
#'          absorb=c("death","death.s"),
#'          cov=c("age"),cov_value=c(0))
#'sim.MTL(obj=fit,state=c(2,3),nsim=1000,L=120)
#'
#'
#'data(SimCml)
#'SimCml=SimCml[,c(1,2,9,10,11,12,13,14)]
#' fit=optim.fit(data=SimCml,
#'          transM=matrix(c(0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0),4,byrow=TRUE),
#'          nstate=4,state_label=c("diagnose","cp1","sct","death"),
#'          event_label=c("cp1.s","sct.s","death.s"),
#'          treatment=c("sct","sct.s"),absorb=c("death","death.s"),
#'          cov=c("age","sex"),cov_value=c(0,1))
#'sim.MTL(obj=fit,nsim=1000,L=120)}
#'
sim.MTL<-function(obj,state=NULL,nsim=1000,L=120)
{
  ## input
  ## obj=obj
  ## path=c(6,7) or path=null(output all paths)
  ## nsim=1000
  ## L=120
  if(is.null(obj$overall)) stop("Please use object returned from optim.fit")

  if(is.null(obj$overall$treatment))
  {
    nstate=obj$overall$nstate
    ntrans=obj$overall$ntrans
    absorbstate=obj$overall$nstate
    tpA=obj$overall$tranProb
    states=obj$overall$state_label


    if(!is.null(state))
    {
      MTL.state=c()
      for(s in state)
      {
        output=sim.los(tpA,s,absorbstate,nstate,ntrans,nsim,L)
        MTL.state=rbind(MTL.state,cbind(state=s,MTL=mean(output$los)))
      }

    }

    if(is.null(state)){
      object=list()
      object
    }else
    {
      object=list()
      object$state.MTL=MTL.state
      object$state.table=data.frame(state=1:obj$overall$nstate,label=obj$overall$state_label)
      object
    }

  }else if(!is.null(obj$overall$treatment))
  {
    ###################################################################################
    ### no treatment
    nstate=obj$no_treat$nstate
    ntrans=obj$no_treat$ntrans
    absorbstate=obj$no_treat$nstate
    tpA=obj$no_treat$tranProb
    states=obj$no_treat$state_label

    MTL.state=c()
    for(s in 1:(nstate-1))
    {
      output=sim.los(tpA,s,absorbstate,nstate,ntrans,nsim,L)
      MTL.state=rbind(MTL.state,cbind(state=s,MTL=mean(output$los)))
    }
    before=list(state.MTL=MTL.state)

    ###  treatment
    nstate=obj$treat$nstate
    ntrans=obj$treat$ntrans
    absorbstate=obj$treat$nstate
    tpA=obj$treat$tranProb
    states=obj$treat$state_label

    MTL.state=c()
    for(s in 1:(nstate-1))
    {
      output=sim.los(tpA,s,absorbstate,nstate,ntrans,nsim,L)
      MTL.state=rbind(MTL.state,cbind(state=s,MTL=mean(output$los)))
    }
    after=list(state.MTL=MTL.state)

    object=list()
    object$strategies=data.frame(straregies=paste0("no treatment at ",obj$no_treat$state_label[-obj$no_treat$nstate]),MTL=before$state.MTL[,2])
    object$strategies=rbind(object$strategies,data.frame(straregies=paste0("treatment at ",obj$treat$state_label[-obj$treat$nstate]),MTL=after$state.MTL[,2]))
    object
  }


}
