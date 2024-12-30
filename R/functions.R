
#' @importFrom stats ecdf
CompKMcurves=function(id,time,status,t){

  # sort the data by id and time
  o=order(id,time)
  id=id[o]
  time=time[o]
  status=status[o]


  #Number of illness states
  K=max(status)-1

  #number of patients
  uid=unique(id)
  n=length(uid)


  m=length(t)


  obj=ecdf(t)
  #n-vector


  Surv=matrix(NA,K+1,m)

  psi=array(NA,c(K+1,n,m))

  for (k in 1:(K+1)){
    tmp=data.frame(id,time,status)[status==0|status>=k,]
    data.comp.k=tmp[!duplicated(tmp[,"id"]),]

    time.rank.t=round(obj(data.comp.k[,"time"])*m)
    delta=(data.comp.k[,"status"]!=0)+0

    ##########################################
    # Same as univariate version from here on
    #########################################

    # Output : n x m matrices
    # Counting and at risk processes
    dN=matrix(0,n,m)
    R=matrix(0,n,m)

    for (i in 1:n){
      index=time.rank.t[i]
      if (index>0)
        R[i,1:index]=1
      if (delta[i]==1){
        dN[i,index]=1
      }
    }

    # dN[1:20,1:20]
    tot.R=colSums(R)
    dL=rep(0,m)


    dL[tot.R>0]=(colSums(dN)/tot.R)[tot.R>0]
    # dL[1:20]

    Surv[k,]=cumprod(1-dL)

    dM=dN-R*matrix(rep(dL,each=n),n,m)

    pi=tot.R/n

    # influence matrix: n x m
    dpsi=matrix(0,n,m)
    mp=sum(pi>0)
    dpsi[,1:mp]=matrix(rep(1/pi[1:mp],each=n),n,mp)*dM[,1:mp]

    psi[k,,]=t(apply(dpsi, 1, cumsum))
  }



  # Surv[,1:100]


  return(list(Surv=Surv,psi=psi,t=t))
}



#' Estimate restricted mean times in favor of treatment
#'
#' Estimate and make inference on the overall and component-wise
#' restricted mean times in favor of treatment.
#'
#' @param id A vector of id variable.
#' @param time A vector of follow-up times.
#' @param status For \code{type="multistate"}, k = entering into state \eqn{k}
#' (\eqn{K+1} represents death) and 0 = censoring;
#' For \code{type="recurrent"}, 1 = recurrent event, 2 = death,
#' and 0 = censoring;
#' @param trt A vector of binary variable for treatment group.
#' @param type \code{"multistate"} = multistate data; \code{"recurrent"} =
#' recurrent event data.
#' @param formula A formula object. For multistate data, use \code{ms(id,time,status)~trt};
#' for recurrent event data, use \code{rec(id,time,status)~trt}.
#' @param data A data frame, which contains the variables names in the formula.
#' @param ... Further arguments.
#' @return An object of class \code{rmtfit}. See \code{\link{rmtfit.object}} for details.
#' @examples
#' #######################
#' # Multistate outcome  #
#' #######################
#' # load the colon cancer trial data
#' library(rmt)
#' head(colon_lev)
#' # fit the data
#' obj=rmtfit(ms(id,time,status)~rx,data=colon_lev)
#' # print the event numbers by group
#' obj
#' # summarize the inference results for tau=7.5 years
#' summary(obj,tau=7.5)
#'
#' ############################
#' # Recurrent event outcome  #
#' ############################
#' # load the HF-ACTION trial data
#' library(rmt)
#' head(hfaction)
#' # fit the data
#' obj=rmtfit(rec(patid,time,status)~trt_ab,data=hfaction)
#' # print the event numbers by group
#' obj
#' # summarize the inference results for tau=3.5 years
#' summary(obj,tau=3.5,Kmax=4) # aggregating results for recurrent-event
#' # frequency >=4.
#'
#'
#' @keywords rmtfit
#' @export
#' @aliases rmtfit
#' @seealso \code{\link{rmtfit.object}},
#' \code{\link{summary.rmtfit}}, \code{\link{plot.rmtfit}}, \code{\link{bouquet}}.
#'
rmtfit=function(...)UseMethod("rmtfit")



#'@describeIn rmtfit Default
#'@export
#'@importFrom stats median
rmtfit.default=function(id,time,status,trt,type="multistate",...){

  if (type=="recurrent"){

    # sort by id and time
    o=order(id,time)

    id=id[o]
    time=time[o]
    status=status[o]
    trt=trt[o]

    freq_rec=table(id)
    id_list=names(freq_rec)
    # recreate the status variable

    K=max(freq_rec)-1


    l_status=status[status!=1]
    l_status[l_status==2]=K+1


    re_status=function(x){
      y=1:freq_rec[x]
      y[freq_rec[x]]=l_status[x]
      return(y)
    }

    status=unlist(lapply(1:length(freq_rec),re_status))
  }





  ###########################################
  # Need a chunk of code for format checking
  ###########################################

  trt.level=sort(unique(trt))

  if (length(trt.level)!=2){
    stop("The explanatory variable (i.e., trt) must be binary!")
  }


  trt0=trt.level[1]
  trt1=trt.level[2]

  # trt=(data$size>0)
  t=sort(unique(time[status>0]))

  ind1=(trt==trt1)
  ind0=(trt==trt0)

  id1=id[ind1]
  time1=time[ind1]
  status1=status[ind1]

  id0=id[ind0]
  time0=time[ind0]
  status0=status[ind0]

  n1=length(unique(id1))
  n0=length(unique(id0))

  n=n1+n0
  K=max(status)-1
  m=length(t)




  ####### Descriptive statistics ######
  # Event frequencies
  desc=table(trt,status)
  # Remove empty categories
  if (any(rowSums(desc)==0)){
    desc=desc[-which(rowSums(desc)==0),]
  }
  # Median follow-up times
  FU1=median(time1[status1==0|status1==K+1])
  FU0=median(time0[status0==0|status0==K+1])

  desc[,1]=c(n0,n1)
  desc=cbind(desc,c(FU0,FU1))
  state=ifelse(type=="multistate","State","Event")

  colnames(desc)=c("N",paste(state,1:K),"Death","Med follow-up time")
  #########Output: desc#######



  # system.time({
  obj1=CompKMcurves(id1,time1,status1,t)
  # })

  # system.time({
  obj0=CompKMcurves(id0,time0,status0,t)
  # })




  # dim(psi1)

  Surv1=obj1$Surv
  psi1=obj1$psi

  Surv0=obj0$Surv
  psi0=obj0$psi

  # dim(Surv0)
  #
  # dim(psi1)
  # dim(psi0)


  Surv1.lead=rbind(Surv1[2:(K+1),],rep(1,m))
  #psi1.diff[k,,]=psi1_{k+1}-psi1_k
  psi1.lead=array(NA,c(K+1,n1,m))
  psi1.lead[1:K,,]=psi1[2:(K+1),,]
  psi1.lead[(K+1),,]=0


  Surv0.lead=rbind(Surv0[2:(K+1),],rep(1,m))
  #psi0.diff[k,,]=psi0_{k+1}-psi0_k
  psi0.lead=array(NA,c(K+1,n0,m))
  psi0.lead[1:K,,]=psi0[2:(K+1),,]
  psi0.lead[(K+1),,]=0


  dt=c(t[1],diff(t))



  ## estimates
  ## Surv.10 and Surv.01 are eta_10 and eta_01
  ## in the paper, respectively
  Surv.10=Surv1*Surv0.lead
  Surv.01=Surv0*Surv1.lead

  Surv.prod=Surv1*Surv0
  eta10=Surv.10-Surv.prod
  eta01=Surv.01-Surv.prod

  eta=Surv.10-Surv.01

  dt.mat=matrix(rep(dt,each=K+1),K+1,m)


  mu10=t(apply(eta10*dt.mat,1,cumsum))
  mu01=t(apply(eta01*dt.mat,1,cumsum))

  # mu10[,1100:1120]
  # mu01[,1100:1120]
  #
  # dim(mu10)

  mu=mu10-mu01

  # mu10[,which.min(t<15)]
  # mu01[,which.min(t<15)]

  nmax=max(n1,n0)

  Surv.10.array=aperm(array(rep(Surv.10,nmax),c(K+1,m,nmax)),c(1,3,2))
  Surv.01.array=aperm(array(rep(Surv.01,nmax),c(K+1,m,nmax)),c(1,3,2))

  psi1.diff=Surv.01.array[,1:n1,]*psi1.lead-Surv.10.array[,1:n1,]*psi1
  psi0.diff=-Surv.10.array[,1:n0,]*psi0.lead+Surv.01.array[,1:n0,]*psi0

  dt.array=array(rep(dt,each=(K+1)*nmax),c(K+1,nmax,m))

  IC1.intgran=psi1.diff*dt.array[,1:n1,]
  IC0.intgran=psi0.diff*dt.array[,1:n0,]

  IC1=apply(IC1.intgran,1:2,cumsum)
  IC0=apply(IC0.intgran,1:2,cumsum)
  # dim(IC1.tot)



  # strangely, the following looping
  # is faster than applying magin sums (commented code)
  IC1.tot=matrix(0,m,n1)
  IC0.tot=matrix(0,m,n0)
  for (k in 1:(K+1)){
    IC1.tot=IC1.tot+IC1[,k,]
    IC0.tot=IC0.tot+IC0[,k,]
  }

  var1=rbind(t(apply(IC1^2,1:2,mean)),rowMeans(IC1.tot^2))/n1
  var0=rbind(t(apply(IC0^2,1:2,mean)),rowMeans(IC0.tot^2))/n0
  # dim(var1)
  #
  # IC1.tot=t(apply(IC1,c(1,3),sum) )#n x m
  # IC0.tot=t(apply(IC0,c(1,3),sum))
  # # dim(IC1.tot)
  # var1=colMeans(IC1.tot^2)/n1
  # var0=colMeans(IC0.tot^2)/n0

  var=var1+var0

  mu=rbind(mu,colSums(mu))

  result=list(t=t,mu=mu,mu10=mu10,mu01=mu01,IC1=IC1,IC0=IC0,var=var,trt1=trt1,trt0=trt0,
              desc=desc,type=type)
  result$call<-match.call()

  class(result)<-"rmtfit"

  return(result)
}



#'@describeIn rmtfit Formula
#'@importFrom stats model.frame
#'@export
rmtfit.formula=function(formula,data,...){
  mf=model.frame(formula=formula,data=data)

  form=mf[,1]

  type=ifelse(class(form)=="rec","recurrent","multistate")

  id=rownames(form)
  # id=data$id


  time=form[,1]
  status=form[,2]
  trt=mf[,2]


  obj=rmtfit.default(id=id,time=time,status=status,trt=trt,type=type,...)
  obj$call<-match.call()
  obj$formula=formula
  return(obj)

}


#'Create a multistate event object
#' @description Create a multistate event object
#' @param id A vector of id variable.
#' @param time A vector of follow-up times.
#' @param status A vector of event type, \code{k} if transitioning to state \eqn{k},
#' 0 if censored and \eqn{K+1} represents death.
#' @return An object of class \code{ms} used as an argument for \code{\link{rmtfit}}.
#'@export
ms=function(id,time,status){

  ss=cbind(time=time,status=status)
  rownames(ss)=id
  class(ss) <- "ms"
  ss
}

#' Create a recurrent event object
#' @description Create a recurrent event object
#' @param id A vector of id variable.
#' @param time A vector of follow-up times.
#' @param status A vector of event type, 1 = recurrent event, 2 = death, and 0 = censoring;
#' @return An object of class \code{rec} used as an argument for \code{\link{rmtfit}}.
#'@export
rec=function(id,time,status){
  ss=cbind(time=time,status=status)
  rownames(ss)=id
  class(ss) <- "rec"
  ss
}






