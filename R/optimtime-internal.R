#' internal
#'
#' internal
#'
#' internal
#'
#' @keywords internal
#'
#'
#'
sim.los<-function(tpA,startstate0=1,absorbstate=5,nstate=5,ntrans=9,n=1000,L=120){
  tpwork=vector("list",nstate)
  for(i in 1:nstate)
  {
    temp=as.matrix(tpA[[i]])
    tpwork[[i]]=cbind(temp[,1],t(apply(temp[,2:(2+nstate-1)],1,cumsum)))
  }
  los.all=rep(0,n)
  path.all=vector("list",n)
  for(i in 1:n)
  {
    startstate=startstate0
    absorbstate=absorbstate
    index0=1
    los=0
    los.path=los
    state.path=startstate
    while(startstate!=absorbstate&los<L)
    {
      u=runif(1)
      ## time point to transit from startstate
      index=which(tpwork[[startstate]][index0:dim(tpwork[[startstate]])[1],startstate+1]<u)[1]+index0-1
      if(is.na(index)){index=dim(tpwork[[startstate]])[1] }
      nextstate=which(tpwork[[startstate]][index,-1]>u)[1]
      if(is.na(nextstate)) break
      los=tpwork[[startstate]][index,1]
      if(los>L) los=L
      startstate=nextstate
      los.path=c(los.path,los)
      state.path=c(state.path,nextstate)
      index0=index
    }
    los.all[i]=los
    path.all[[i]]=state.path
  }
  return(list(los=los.all,path=path.all))
}
