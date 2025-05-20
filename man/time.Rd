\name{Time functions}
\alias{nanoTime}
\title{
Time functions
}

\description{
Functions to get the current time in different units.
}

\usage{
nanoTime()
}

\arguments{
}

\details{
    
\itemize{
\item \code{nanoTime}: The current nanoseconds of the system.
}
}


\value{
Each function returns a numeric value.
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{runif}}, \code{\link{rbeta}}, \code{\link{rexp}}, \code{\link{rchisq}}, \code{\link{rgamma}}, \code{\link{rcauchy}}, \code{\link{rt}}, \code{\link{rgeom}}
}

\examples{
x <- nanoTime()
}
