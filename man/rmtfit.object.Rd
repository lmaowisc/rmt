\name{rmtfit.object}
\alias{rmtfit.object}

\title{
  Estimated restricted mean times in favor of treatment
}
\description{
This class of objects is returned by the \code{rmtfit} class of functions.
Objects of this class have methods for the functions \code{print},
\code{summary}, \code{plot}, and \code{bouquet}.
}
\value{
\item{t}{A vector of follow-up times \eqn{\tau}.
}
\item{mu}{A matrix with \eqn{K+2} rows; The \eqn{k}th row \eqn{(k=1,\ldots, K)}
is \eqn{\mu_k(\tau)}, the restricted mean time in favor of treatment on the
\eqn{k}th state (or recurrent event); The \eqn{(K+1)}th row is the
net restricted mean survival time; The last row is the overal effect \eqn{\mu(\tau)}.}

\item{var}{A matrix with \eqn{K+2} rows containing the variance estimates
for \code{mu}.
}

\item{mu10, mu01}{Matrices with \eqn{K+1} rows; The \eqn{k}th row \eqn{(k=1,\ldots, K)}
is the restricted mean win (\code{mu10}) and loss (\code{mu01}) times by the treatment on the
\eqn{k}th state (or recurrent event); The \eqn{(K+1)}th row is that
on the survival time.}

\item{...}{}
}

\seealso{
\code{\link{rmtfit}}, \code{\link{summary.rmtfit}}, \code{\link{plot.rmtfit}}, \code{\link{bouquet}}.
}
\keyword{rmtfit}
