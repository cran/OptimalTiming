## optim.fit
#' Fit multi-state model for optimization
#'
#' This function produces transition probabilities for given covariates values in multi-state models.
#'
#' For optim.fit, transition probabilities are estimated under Markov assumption, which implies that
#' the probability of transition to a future state depends only on the present state, not on the history.
#' Taking covariates at baseline into consideration, transition probabilities can be subject-specific.
#' Cox proportional hazards model is used to fit transition hazards among multiple states by assuming
#' each transition has its own baseline hazard, and covariates have different effects on different transitions.
#'
#' Let \eqn{\mathbf{S}={1,2,\cdots,S}} denote the states in the multi-state model
#'  and \eqn{X(t)} be a random process taking values from \eqn{\mathbf{S}}.
#' Denote \eqn{\alpha_{gh}(t)} as hazard ratio or transition intensity and \eqn{Z} as baseline covariates.
#' The instantaneous risk of a transition from state \eqn{g} into state \eqn{h} at time \eqn{t} can be
#' fitted by semi-parametric Cox model:
#' \deqn{\alpha_{gh}(t|Z)=\alpha_{gh,0}exp(\beta^{T}Z_{gh}).}
#' The cumulative hazard ratio is defined as \eqn{A_{gh}(t)=\int_0^t \alpha_{gh}(u)du}.
#' Primary interest in this function is to estimate transition probability \eqn{P_{gh}(s,t)=P(X(t)=h|X(s)=g)},
#' indicating the chance of transition from state \eqn{g} at time \eqn{s} to state \eqn{h} at time \eqn{t}.
#' Written in matrix form, transition probability matrix \eqn{\mathbf{P}(t)} can be calculated by
#' means of a product integral: \eqn{\mathbf{P}(s,t)=\prod_{(s,t]}(\mathbf{I}+d\mathbf{A}(u))},
#' where \eqn{\mathbf{A}(t)} is a transition intensity matrix. Both \eqn{\mathbf{P}} and \eqn{\mathbf{A}}
#' are \eqn{S \times S} matrix.
#'
#' The data format required by this function is wide format, which can be regarded as the augmented
#' data used in single event survival analysis. For example, if there is a "recurrence" state in a
#' multi-state model, two variable are needed to describe this event, namely, "rec" and "rec.s".
#' The former is a time variable, indicating the time from initiation of the study
#' to the occurrence of this state, while the latter is an indicator variable with 1 for occurrence and 0 for
#' censoring.
#' If the event is censored for some patients, use the maximum follow-up instead of the event time.
#' Other states are prepared in the same way. Thus, each row in the augmented data
#' summarize all possible events for a single subject. For covariates, they are located at the end
#' of each row.
#'
#'
#' If the time of new treatment initiation is provided in data, the argument \code{treatment}
#'  should be assigned as, eg.\code{treatment=c("sct","sct.s")}. Additionally, the argument
#'  \code{state_lable} and \code{event_label} shoud be arranged in such order: pre-treatment state,
#'  treatment state, post-treatment states and  absorbing state. Assume treatment may take place
#'  at any pre-treatment states.  In this case,
#'  \code{optim.fit} function automatically fit two multistate models, one for post-treatment states
#'   if a new treatment is carried out , and
#' the other  for pre-treatment states if a new treatment is not carried out.
#' Thus, comparison among strategies of whether and when to initiate new treatment can be
#' performed in \code{sim.MTL} function.   If \code{treatment=NULL},
#' a single multistate model will be fitted.
#'
#'
#'
#' @param data Data frame in wide format where each row in the data corresponds to a single subject.
#'  Time to a state and the occurrence of the state come in pairs. If a state is not occur,
#'  use the time to an absorbing state or censoring time instead. Covariates are added at the end
#'  of each row.
#' @param transM A \eqn{nstate \times nstate} matrix used to indicate the transitions in multi-state model.
#' If a transition exists between two states,
#' set 1 in a corresponding location, otherwise set 0.
#' @param nstate Number of states incorporated in the multi-state model.
#' @param state_label A character vector of length \code{nstate}  containing the names of states.
#' The elements in \code{state_label} are
#' extracted from the column names of \code{data}, except for the first one, which is a potential state
#' at the initiation of a study for each subject. Assume all subjects have the same initial state.
#' @param event_label A character vector of length \code{nstate-1}, indicating the occurrence of each state.
#' The first state in \code{state_label} do not need an indiator, as it always exists.
#' @param treatment A character vector of length 2, indicating whether there is a treatment variable
#' available. If true, the name and indicator of this treatment extracted from \code{state_label}
#' and \code{event_label} consist of \code{treatment}. If not, \code{treatment}=NULL. See details.
#' The default value is NULL.
#' @param absorb A character vector of length 2, indicating the name and indicator of the absorb state.
#' @param cov A character vector containing the names of covariates that have some effect to transition
#' probabilities.
#' @param cov_value A numeric vector containing the values of covariates corresponding to \code{cov}.
#' \code{cov_value} are used to calculated subject specified transition probabilities.
#' @return If \code{treatment} is NULL, a list object called "overall" is output with elements:
#' \item{transMat}{A transition matrix describing the states and transitions in multi-state model.}
#' \item{tranProb}{A list of size \code{nstate} recording the transition probabilities form each state to  another along with standard errors. Element \code{[[s]]} is a data frame containing transition  probabilities from state \code{s} to state \code{1,2,...,nstate}.These transition probabilities are  time-varying over distinct transition time points.}
#' \item{coxobj}{An object returned by \code{coxph()} function in \code{survival} package.}
#' \item{ntrans}{The number of available transitions among multiple states.}
#' \item{...}{Other variables that extracted from the original input.}
#' If \code{treatment} is not NULL, three list objects called "overall", "treat","no_treat" are output.
#' A list "overall" contains elements:
#' \item{transMat}{A transition matrix describing the states and transitions in multi-state model.}
#' \item{ntrans}{The number of available transitions among multiple states.}
#' \item{...}{Other variables that extracted from the original input.}
#' A list "no_treat" contains elements by fitting a model for pre-treatment states:
#' \item{transMat}{A transition matrix describing the states and transitions if the new treatment is  not carried out.}
#' \item{tranProb}{A list recording the transition probabilities among pre-treatment states.}
#' \item{coxobj}{An object returned by \code{coxph()} function in \code{survival} package.}
#' \item{ntrans}{The number of available transitions among pre-treatment states.}
#' \item{data}{A data set contaning states if the new treatment is not carried out.}
#' \item{nstate}{The number of pre-treatment states.}
#' \item{...}{Other variables that extracted from the original input.}
#' A list "treat" contains elements by fitting a model if a new treatment is carried out:
#' \item{transMat}{A transition matrix describing the states and transitions if the new treatment is  carried out.}
#' \item{tranProb}{A list recording the transition probabilities among post-treatment states.}
#' \item{coxobj}{An object returned by \code{coxph()} function in \code{survival} package.}
#' \item{ntrans}{The number of available transitions among post-treatment states.}
#' \item{data}{A data set contaning states if the new treatment is carried out.}
#' \item{nstate}{The number of post-treatment states.}
#' \item{...}{Other variables that extracted from the original input.}
#'
#'
#'
#' @aliases optim.fit
#' @references de Wreede LC, Fiocco M, and Putter H (2010). The mstate package for estimation and prediction in non- and semi-parametric multi-state and competing risks models. Computer Methods and Programs in Biomedicine 99, 261–274.
#' @export
#' @examples
#'
#'\dontrun{
#' library(OptimalTiming)
#' ## read data
#' data(SimCml)
#'
#' ## fit multistate model if the time to new treatment initiation is available in SimCml
#' fit=optim.fit(data=SimCml,
#'          transM=matrix(c(0,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,1,1,1,
#'          0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0),7,byrow=TRUE),
#'          nstate=7,state_label=c("diagnose","cp1","ap","cp2","bc","sct","death"),
#'          event_label=c("cp1.s","ap.s","cp2.s","bc.s","sct.s","death.s"),
#'          treatment=c("sct","sct.s"),absorb=c("death","death.s"),
#'          cov=c("age"),cov_value=c(0))
#'
#' ## view the content of this object
#' objects(fit)
#'
#' ## output transition probabilities
#' fit$treat$tranProb
#' fit$no_treat$tranProb}
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
#'
#'
#'data(SimCml)
#'SimCml=SimCml[,c(1,2,9,10,11,12,13,14)]
#' fit=optim.fit(data=SimCml,
#'          transM=matrix(c(0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0),4,byrow=TRUE),
#'          nstate=4,state_label=c("diagnose","cp1","sct","death"),
#'          event_label=c("cp1.s","sct.s","death.s"),
#'          treatment=c("sct","sct.s"),absorb=c("death","death.s"),
#'          cov=c("age","sex"),cov_value=c(0,1))}

optim.fit <- function(data,transM,nstate,state_label,event_label,treatment=NULL,absorb,cov,cov_value){
  ## input
  ## data:framework
  ## transM=matrix(c(0,1,1,1,1,1,1,0,0,1,0,1,1,1,0,0,0,1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0),nstate,byrow=T)
  ## nstate:nstate=7
  ## state_label:state_label=c("diagnose","cp1","ap","cp2","bc","sct","death")
  ## event_label:event_label=c("cp1.s","ap.s","cp2.s","bc.s","sct.s","death.s")
  ## treatment=c("sct","sct.s")
  ## absorb=c("death","death.s")
  ## cov:cov=c("agegp","sex")
  ## cov_value:cov_value=c(0,1)

  if(is.null(treatment))
  {
    transM_temp=vector("list",nstate)
    for(i in 1:nstate){ transM_temp[[i]]=which(transM[i,]==1)}

    ## transition matrix
    tmat <- mstate::transMat(x=transM_temp,names=state_label)

    N.trans=max(tmat,na.rm=T)

    ## arrange data in long format
    datalong <- mstate::msprep(time=c(NA,state_label[-1]),status=c(NA,event_label),data=data,keep=names(data)[(2*nstate-1):dim(data)[2]],trans=tmat)

    ## expanded covariates
    datalongcov <- mstate::expand.covs(datalong,cov)

    ## Cox model with different covariate
    cov.reg=c()
    for(i in 1:length(cov)){
      cov.reg=paste0(cov.reg,"+",paste0(cov[i],".",1:N.trans,collapse='+'))}

    formula=paste("Surv(Tstart,Tstop,status)","~", substr(cov.reg,2,10000000),"+","strata(trans)",sep=' ')

    cx <- coxph(as.formula(formula),data=datalongcov,method="breslow")

    ## estimate transition intensity for new patient
    covcol=which(colnames(datalong)%in%cov)
    condition=which(datalong[,covcol]==cov_value)[1]
    newdata=datalong[rep(condition,N.trans),]
    newdata$trans=1:N.trans
    attr(newdata,"trans")<-tmat
    newdata=mstate::expand.covs(newdata,cov)
    newdata$strata=newdata$trans
    msfA=mstate::msfit(cx,newdata,trans=tmat)

    ## estimate transition probability for new patient
    tpA=mstate::probtrans(msfA,predt=0,direction="fo")

    ##### output
    object=list()
    object$overall$transMat=tmat
    object$overall$tranProb=tpA[1:nstate]
    object$overall$coxobj=cx
    object$overall$ntrans=N.trans
    object$overall$data=data
    object$overall$transM=transM
    object$overall$nstate=nstate
    object$overall$state_label=state_label
    object$overall$event_label=event_label
    object$overall$treatment=treatment
    object$overall$absorb=absorb
    object$overall$cov=cov
    object$overall$cov_value=cov_value
    return(object)

  }else if(!is.null(treatment))
  {
    ##overall
    transM_temp=vector("list",nstate)
    for(i in 1:nstate){ transM_temp[[i]]=which(transM[i,]==1)}

    ## transition matrix
    tmat <- mstate::transMat(x=transM_temp,names=state_label)

    N.trans=max(tmat,na.rm=T)

    ## before transplant
    ind_before=1:(which(names(data)==treatment[1])-1)
    ind_before=c(ind_before,which( names(data)%in% absorb))
    #ind_before=c(ind_before,which( names(data)%in% cov))
    ind_before=c(ind_before,(which(names(data)==absorb[1])+2):dim(data)[2])
    data1=data[data[,which(names(data)==treatment[2])]==0,ind_before]

    ind_before_tranM=1:(which(state_label==treatment[1])-1)
    ind_before_tranM=c(ind_before_tranM,which(state_label==absorb[1]))
    transM1=transM[ind_before_tranM,ind_before_tranM]
    nstate1=dim(transM1)[1]

    state_label1=state_label[ind_before_tranM]
    event_label1=event_label[ind_before_tranM-1]

    transM_temp1=vector("list",nstate1)
    for(i in 1:nstate1){ transM_temp1[[i]]=which(transM1[i,]==1)}

    ## transition matrix
    tmat1 <- mstate::transMat(x=transM_temp1,names=state_label1)

    N.trans1=max(tmat1,na.rm=T)

    ## arrange data in long format
    datalong1 <- mstate::msprep(time=c(NA,state_label1[-1]),status=c(NA,event_label1),data=data1,keep=names(data1)[which(names(data)==absorb[1]):dim(data1)[2]],trans=tmat1)

    ## expanded covariates
    datalongcov1 <- mstate::expand.covs(datalong1,cov)

    ## Cox model with different covariate
    cov.reg1=c()
    for(i in 1:length(cov)){
      cov.reg1=paste0(cov.reg1,"+",paste0(cov[i],".",1:N.trans1,collapse='+'))}

    formula1=paste("Surv(Tstart,Tstop,status)","~", substr(cov.reg1,2,10000000),"+","strata(trans)",sep=' ')

    cx1 <- coxph(as.formula(formula1),data=datalongcov1,method="breslow")

    ## estimate transition intensity for new patient
    covcol1=which(colnames(datalong1)%in%cov)
    condition1=which(datalong1[,covcol1]==cov_value)[1]
    newdata1=datalong1[rep(condition1,N.trans1),]
    newdata1$trans=1:N.trans1
    attr(newdata1,"trans")<-tmat1
    newdata1=mstate::expand.covs(newdata1,cov)
    newdata1$strata=newdata1$trans
    msfA1=mstate::msfit(cx1,newdata1,trans=tmat1)

    ## estimate transition probability for new patient
    tpA1=mstate::probtrans(msfA1,predt=0,direction="fo")


    ##############################################################
    ## after transplant
    ind_treat=which(names(data)==treatment[1])
    ind_death=which(names(data)==absorb[1])
    ind_after=1:(which(names(data)==treatment[1])-1)
    data2=data[data[,names(data)==treatment[2]]==1,]
    for(i in 1:dim(data2)[1])
    {
      data2[i,which(data2[i,ind_after]==1)-1]=data2[i,ind_treat]
      data2[i,which(data2[i,ind_after]==0)-1]=data2[i,ind_death]
    }
    ind_after=c(ind_after,(ind_treat+2):(ind_death+1))
    ind_after=c(ind_after,which(names(data2)%in%cov))
    data2=data2[,ind_after]


    ind_after_tranM=1:(which(state_label==treatment[1])-1)
    nstate2=length(state_label)-1
    transM2=matrix(NA,nstate2,nstate2)
    transM2[1,ind_after_tranM]=c(0,rep(1,length(ind_after_tranM)-1))
    transM2[ind_after_tranM[-1],ind_after_tranM]=0

    ind_treat_tranM=which(state_label==treatment[1])
    ind_death_tranM=which(state_label==absorb[1])
    if(ind_treat_tranM+1==ind_death_tranM)
    {
      transM2[,nstate2]=transM[-ind_treat_tranM,-ind_treat_tranM][,nstate2]
      transM2[nstate2,]=transM[-ind_treat_tranM,-ind_treat_tranM][nstate2,]
    }else
    {
      ind_rec.tranM=(ind_treat_tranM+1):(ind_death_tranM-1)
      transM2[1,ind_rec.tranM-1]=0
      transM2[ind_after_tranM[-1],ind_rec.tranM-1]=matrix(rep(transM[ind_treat_tranM,ind_rec.tranM],length(ind_after_tranM)-1),nrow=length(ind_after_tranM)-1,byrow=T)
      transM2[ind_rec.tranM-1,ind_after_tranM]=0
      transM2[ind_rec.tranM-1,ind_rec.tranM-1]=transM[ind_rec.tranM,ind_rec.tranM]
      transM2[,nstate2]=transM[-ind_treat_tranM,-ind_treat_tranM][,nstate2]
      transM2[nstate2,]=transM[-ind_treat_tranM,-ind_treat_tranM][nstate2,]
    }

    state_label2=state_label[-ind_treat_tranM]
    event_label2=event_label[-(ind_treat_tranM-1)]

    transM_temp2=vector("list",nstate2)
    for(i in 1:nstate2){ transM_temp2[[i]]=which(transM2[i,]==1)}

    ## transition matrix
    tmat2 <- mstate::transMat(x=transM_temp2,names=state_label2)

    N.trans2=max(tmat2,na.rm=T)

    ## arrange data in long format
    datalong2 <- mstate::msprep(time=c(NA,state_label2[-1]),status=c(NA,event_label2),data=data2,keep=names(data2)[which(names(data)==absorb[1]):dim(data2)[2]],trans=tmat2)

    ## expanded covariates
    datalongcov2 <- mstate::expand.covs(datalong2,cov)

    ## Cox model with different covariate
    cov.reg2=c()
    for(i in 1:length(cov)){
      cov.reg2=paste0(cov.reg2,"+",paste0(cov[i],".",1:N.trans2,collapse='+'))}

    formula2=paste("Surv(Tstart,Tstop,status)","~", substr(cov.reg2,2,10000000),"+","strata(trans)",sep=' ')

    cx2 <- coxph(as.formula(formula2),data=datalongcov2,method="breslow")

    ## estimate transition intensity for new patient
    covcol2=which(colnames(datalong2)%in%cov)
    condition2=which(datalong2[,covcol2]==cov_value)[1]
    newdata2=datalong2[rep(condition2,N.trans2),]
    newdata2$trans=1:N.trans2
    attr(newdata2,"trans")<-tmat2
    newdata2=mstate::expand.covs(newdata2,cov)
    newdata2$strata=newdata2$trans
    msfA2=mstate::msfit(cx2,newdata2,trans=tmat2)

    ## estimate transition probability for new patient
    tpA2=mstate::probtrans(msfA2,predt=0,direction="fo")



    ##### output
    object=list()
    object$no_treat$transMat=tmat1
    object$no_treat$tranProb=tpA1[1:nstate1]
    object$no_treat$coxobj=cx1
    object$no_treat$ntrans=N.trans1
    object$no_treat$data=data1
    object$no_treat$transM=transM1
    object$no_treat$nstate=nstate1
    object$no_treat$state_label=state_label1
    object$no_treat$event_label=event_label1
    object$no_treat$cov=cov
    object$no_treat$cov_value=cov_value

    object$treat$transMat=tmat2
    object$treat$tranProb=tpA2[1:nstate2]
    object$treat$coxobj=cx2
    object$treat$ntrans=N.trans2
    object$treat$data=data2
    object$treat$transM=transM2
    object$treat$nstate=nstate2
    object$treat$state_label=state_label2
    object$treat$event_label=event_label2
    object$treat$cov=cov
    object$treat$cov_value=cov_value

    object$overall$transMat=tmat
    object$overall$ntrans=N.trans
    object$overall$data=data
    object$overall$transM=transM
    object$overall$nstate=nstate
    object$overall$state_label=state_label
    object$overall$event_label=event_label
    object$overall$treatment=treatment
    object$overall$absorb=absorb
    object$overall$cov=cov
    object$overall$cov_value=cov_value

    return(object)

  }

}
