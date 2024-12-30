
#' Print a short summary of rmtfit objects
#'
#' Print the results for the restricted mean times in favor of treatment.
#'
#' @param x An object returned by \code{\link{rmtfit}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#' @export
print.rmtfit=function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$desc)
}


#' Summary of the analysis results
#'
#' Summarize the overall and stage-wise inferential results for the restricted
#' mean times in favor of treatment at a user-specified length of follow-up.
#'
#'
#' @param object An object returned by \code{\link{rmtfit}}.
#' @param tau A positive real number for the follow-up time; Default is the maximum
#' event time in the data.
#' @param Kmax A positive integer; If specified, the stage-wise estimates over
#' \code{Kmax}\eqn{,\ldots,K} will be aggregated.
#' @param ...  Additional arguments affecting the summary produced.
#'
#' @return An object of class \code{summary.rmtfit} with components
#' \item{WL}{A \eqn{2\times(K+2)}-dimensional matrix; Each row contains the estimates for the
#' stage-wise and overall restricted mean win times for each group.}
#' \item{tab}{A \eqn{(K+2)\times 4}-dimensional matrix summarizing the inferential results
#' for the stage-wise and overall restricted mean times in favor of treatment; Columns include
#' \code{Estimate}, \code{Std.Err}, \code{Z value}, and \code{Pr(>|z|)}.}
#' \item{...}{}
#' @seealso \code{\link{rmtfit}}, \code{\link{plot.rmtfit}}, \code{\link{bouquet}}.
#' @keywords rmtfit
#' @importFrom stats pchisq
#' @examples
#' #See examples for rmtfit().
#' @export
summary.rmtfit=function(object,tau=NULL,Kmax=NULL,...){

  x <- object

  t=x$t
  m=length(t)
  t.min=t[1]
  t.max=t[m]

  #### check on the input value tau ##########
  # If null, set tau=t.max
  if (is.null(tau)){
    tau=t[m]
  }else{
    if (tau<t.min||tau>t.max){
      stop(
        paste0("tau value is outside the range of empirical observations.\n"),
        "Specify a tau within [",t.min,", ",t.max,"].\n")
    }
  }
  ###########################################

  # index of tau in t
  m0=sum(t<=tau)

  mu.tau=x$mu[,m0]
  K=length(mu.tau)-2


  mu10.tau=x$mu10[,m0]
  mu01.tau=x$mu01[,m0]
  var.tau=x$var[,m0]
  trt1=as.character(x$trt1)
  trt0=as.character(x$trt0)

  # winning and losing time on each state
  WL=matrix(NA,2,K+2)
  rownames(WL)=c(trt0,trt1)
  state=ifelse(x$type=="multistate","State","Event")

  colnames(WL)=c(paste(state,1:K),"Survival","Overall")

  WL[1,]=c(mu01.tau,sum(mu01.tau))
  WL[2,]=c(mu10.tau,sum(mu10.tau))

  # Inference table on the total and component-wise RMT



  if (is.null(Kmax)){
    est=mu.tau
    se=sqrt(var.tau)
    pval=1-pchisq((est/se)^2,1)
    main.tab=cbind(Estimate = est,
                   StdErr = se,
                   z.value = est/se,
                   p.value = pval)
    rownames(main.tab)=c(paste(state,1:K),"Survival","Overall")
  }else{
    if (Kmax<1||Kmax>K){
      stop(paste0("Kmax must be between 1 and ", K,"!"))
    }else{
      if (Kmax>1){
        est=c(mu.tau[1:(Kmax-1)],sum(mu.tau[Kmax:K]),mu.tau[(K+1):(K+2)])
        se=sqrt(c(var.tau[1:(Kmax-1)],sum(var.tau[Kmax:K]),var.tau[(K+1):(K+2)]))
      }else{
        est=c(sum(mu.tau[1:K]),mu.tau[(K+1):(K+2)])
        se=sqrt(c(sum(var.tau[1:K]),var.tau[(K+1):(K+2)]))
      }

      pval=1-pchisq((est/se)^2,1)
      main.tab=cbind(Estimate = est,
                     StdErr = se,
                     z.value = est/se,
                     p.value = pval)
      if (Kmax==1){
        rownames(main.tab)=c(paste(state,"1+"),"Survival","Overall")
      }else{
        rownames(main.tab)=c(paste(state,c(1:(Kmax-1),paste0(Kmax,"+"))),"Survival","Overall")
      }
    }
  }


  colnames(main.tab)=c("Estimate", "Std.Err", "Z value", "Pr(>|z|)")



  result=list(WL=WL,tab=main.tab,tau=tau,call=x$call,type=x$type)

  class(result)<-"summary.rmtfit"
  return(result)
}



#' Print method for summary.rmtfit objects
#'
#' Produces a printed summary of the results for the restricted
#' mean times in favor of treatment
#'
#' @param x An object returned by \code{\link{summary.rmtfit}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#'@export
#'@importFrom stats printCoefmat
print.summary.rmtfit=function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n")
  tau=x$tau
  type=x$type


  WL=x$WL
  trt1=rownames(WL)[2]
  cat("Restricted mean winning time by tau = ", round(tau,4), ":\n",sep="")
  print(WL)
  cat("\n")

  tab=x$tab
  cat("Restricted mean time in favor of group ")
  cat(paste0("\"",trt1, "\" by time tau = ",round(tau,4), ":\n"))

  printCoefmat(tab, P.values=TRUE, has.Pvalue=TRUE)

}


#' Plot the estimated treatment effect curve
#'
#' Plot the estimated overall or stage-wise restricted mean times in favor of treatment as a
#' function of follow-up time.
#' @param x An object returned by \code{\link{rmtfit}}.
#' @param k If specified, \eqn{\mu_k(\tau)} is plotted; otherwise, \eqn{\mu(\tau)} is plotted.
#' @param conf If TRUE, 95\% confidence limits for the target curve are overlaid.
#' @param main A main title for the plot
#' @param xlim The x limits of the plot.
#' @param ylim The y limits of the plot.
#' @param xlab A label for the x axis, defaults to a description of x.
#' @param ylab A label for the y axis, defaults to a description of y.
#' @param conf.col Color for the confidence limits if \code{conf=TRUE}.
#' @param conf.lty Line type for the confidence limits if \code{conf=TRUE}.
#' @param ... Other arguments that can be passed to the underlying \code{plot} method.
#' @return No return value, called for side effects.
#' @seealso \code{\link{rmtfit}}, \code{\link{summary.rmtfit}}, \code{\link{bouquet}}.
#' @keywords rmtfit
#' @importFrom stats stepfun qnorm
#' @importFrom graphics abline lines
#' @export
#' @examples
#' # load the colon cancer trial data
#' library(rmt)
#' head(colon_lev)
#' # fit the data
#' obj=rmtfit(ms(id,time,status)~rx,data=colon_lev)
#' # plot overal effect mu(tau)
#' plot(obj)
#' # set-up plot parameters
#' oldpar <- par(mfrow = par("mfrow"))
#' par(mfrow=c(1,2))
#' # Plot of component-wise RMT in favor of treatment over time
#' plot(obj,k=2,conf=TRUE,col='red',conf.col='blue', xlab="Follow-up time (years)",
#'     ylab="RMT in favor of treatment (years)",main="Survival")
#' plot(obj,k=1,conf=TRUE,col='red',conf.col='blue', xlab="Follow-up time (years)",
#'     ylab="RMT in favor of treatment (years)",main="Pre-relapse")
#' par(oldpar)
plot.rmtfit=function(x,k=NULL,conf=FALSE,main=NULL,xlim=NULL, ylim=NULL,xlab="Follow-up time",
  ylab="Restricted mean time in favor",conf.col="black",conf.lty=3,...){
  t=x$t
  mu=x$mu
  K=nrow(mu)-2
  var=x$var

  if (is.null(k)){
    muk=mu[K+2,]
    sek=sqrt(var[K+2,])
  }else{
    if (k%%1==0&&k>=1&&k<=(K+1)){
      muk=mu[k,]
      sek=sqrt(var[k,])
    }else{
      stop(paste0("k must be an integer from 1 to ", K+1,".
  (The default plot is for the total restricted time in favor.)"))
    }
  }


  muk.step=stepfun(t,c(0,muk))


  if (is.null(main)){
    main=paste(x$trt1,"versus",x$trt0)
  }


  if (!conf){
    plot(muk.step,do.points=F,xlab=xlab,ylab=ylab,main=main,...)
    abline(h=0)
  }else{
    za=qnorm(0.975)
    muk.up=muk+za*sek
    muk.lo=muk-za*sek

    if (is.null(ylim)){
      ylim=c(min(muk.lo),max(muk.up))
    }
    plot(muk.step,do.points=F,xlab=xlab,ylab=ylab,main=main,ylim=ylim,...)
    abline(h=0)
    lines(stepfun(t,c(0,muk.up)),col=conf.col,lty=conf.lty,do.points=FALSE)
    lines(stepfun(t,c(0,muk.lo)),col=conf.col,lty=conf.lty,do.points=FALSE)
    # lines(t,muk.up,col=conf.col,lty=conf.lty,do.points=F,lwd=2)
    # lines(t,muk.lo,col=conf.col,lty=conf.lty,do.points=F,lwd=2)
  }
}






#' Bouquet plot
#'
#' Construct the bouquet plot based on the estimated stage-wise restricted mean win/loss
#' times.
#'
#' @param x An object returned by \code{\link{rmtfit}}.
#' @param Kmax A positive integer; If specified, the stage-wise estimates over
#' \code{Kmax}\eqn{,\ldots,K} will be aggregated.
#' @param xlim The x limits of the plot.
#' @param ylim The y limits of the plot.
#' @param xlab A label for the x axis, defaults to a description of x.
#' @param ylab A label for the y axis, defaults to a description of x.
#' @param group.label If \code{TRUE}, group labels will appear on the
#' two sides of the plot.
#' @param cex.group Font size of the group labels if \code{group.label=TRUE}.
#' @param ... Other arguments that can be passed to the underlying \code{plot} method.
#' @return No return value, called for side effects.
#' @seealso  \code{\link{rmtfit}}, \code{\link{summary.rmtfit}},
#'  \code{\link{plot.rmtfit}}.
#' @importFrom grDevices gray
#' @importFrom graphics abline lines polygon text
#' @export
#' @keywords rmtfit
#' @examples
#' # load the colon cancer trial data
#' library(rmt)
#' head(colon_lev)
#' # fit the data
#' obj=rmtfit(ms(id,time,status)~rx,data=colon_lev)
#' # bouquet plot
#' bouquet(obj)
bouquet=function(x,Kmax=NULL,xlim=NULL,ylim=NULL,
                 xlab="Restricted mean win/loss time",ylab="Follow-up time",
                 group.label=TRUE,cex.group=1,...){

  t=x$t
  mu10=x$mu10
  mu01=x$mu01
  K=nrow(mu10)-1



  if (!is.null(Kmax)){
    if (Kmax<1||Kmax>K){
      stop(paste0("Kmax must be between 1 and ", K,"!"))
    }else{
      if (Kmax==1){
        mu10=rbind(colSums(mu10[1:K,]),mu10[(K+1),])
        mu01=rbind(colSums(mu01[1:K,]),mu01[(K+1),])
      }else{
        mu10=rbind(mu10[1:(Kmax-1),],colSums(mu10[Kmax:K,]),mu10[K+1,])
        mu01=rbind(mu01[1:(Kmax-1),],colSums(mu01[Kmax:K,]),mu01[K+1,])
      }
      K=Kmax
    }

  }



  m=length(t)

  xmax=max(colSums(mu10),colSums(mu01))
  if (is.null(xlim)){
    xlim=c(-xmax,xmax)
  }
  ymax=max(t)

  if (group.label){
    y.height=1.05*ymax
  }else{
    y.height=ymax
  }
  if (is.null(ylim)){
    ylim=c(0,y.height)
  }


  plot(rep(0,m),t,lty=1,type='l',ylab=ylab,xlim=xlim,ylim=ylim,
       xlab=xlab,...)

  left.bound=-mu01[K+1,]
  right.bound=mu10[K+1,]
  polygon(c(left.bound,rev(right.bound)),c(t,rev(t)),col=gray(level=0.1),border=NA)

  for (j in 1:K){


    polygon(c(left.bound,rev(left.bound-mu01[K+1-j,])),c(t,rev(t)),col=gray(level=0.1+j*0.6/K),border=NA)
    left.bound=left.bound-mu01[K+1-j,]

    polygon(c(right.bound,rev(right.bound+mu10[K+1-j,])),c(t,rev(t)),col=gray(level=0.1+j*0.6/K),border=NA)
    right.bound=right.bound+mu10[K+1-j,]

  }



  abline(v=0,lwd=2)

  if (group.label){
    text(-0.5*xmax,1.04*ymax,paste("Group",x$trt0),cex=cex.group)
    text(0.5*xmax,1.04*ymax,paste("Group",x$trt1),cex=cex.group)
  }



}



