\name{Random values simulation}
\alias{Runif}
\alias{Rbeta}
\alias{Rexp}
\alias{Rchisq}
\alias{Rgamma}
\alias{Rcauchy}
\alias{Rt}
\alias{Rgeom}
\alias{Rpareto}
\alias{Rfrechet}
\alias{Rlaplace}
\alias{Rgumble}
\alias{Rarcsine}
\alias{Rnorm}
\alias{Runif.mat}
\alias{Rbeta.mat}
\alias{Rexp.mat}
\alias{Rchisq.mat}
\alias{Rgamma.mat}
\alias{Rcauchy.mat}
\alias{Rt.mat}
\alias{Rgeom.mat}
\alias{Rpareto.mat}
\alias{Rfrechet.mat}
\alias{Rlaplace.mat}
\alias{Rgumble.mat}
\alias{Rarcsine.mat}
\alias{Rnorm.mat}
\alias{colRunif}
\alias{colRbeta}
\alias{colRexp}
\alias{colRchisq}
\alias{colRgamma}
\alias{colRgeom}
\alias{colRcauchy}
\alias{colRt}
\alias{colRpareto}
\alias{colRfrechet}
\alias{colRlaplace}
\alias{colRgumble}
\alias{colRarcsine}
\alias{colRnorm}
\alias{set_seed}
\title{
Random values simulation from various distributions
}

\description{
Functions to simulate random values from different probability distributions: uniform, beta, exponential, chi-squared, gamma, Cauchy, t-distribution, and geometric.
}

\usage{
Runif(n, min = 0, max = 1) 
Rbeta(n, alpha, beta)
Rexp(n, rate = 1)
Rchisq(n, df)
Rgamma(n, shape, rate = 1)
Rcauchy(n, location = 0, scale = 1)
Rt(n, df, ncp)
Rgeom(n, prob)
Rpareto(n, shape = 1, scale = 1)
Rfrechet(n, lambda = 1, mu = 0, sigma = 1)
Rlaplace(n, mu = 0, sigma = 1)
Rgumble(n, mu = 0, sigma = 1)
Rarcsine(n, min = 0, max = 1)
Rnorm(n)
Runif.mat(nrow, ncol, min = 0, max = 1) 
Rbeta.mat(nrow, ncol, alpha, beta)
Rexp.mat(nrow, ncol, rate = 1)
Rchisq.mat(nrow, ncol, df)
Rgamma.mat(nrow, ncol, shape, rate = 1)
Rcauchy.mat(nrow, ncol, location = 0, scale = 1)
Rt.mat(nrow, ncol, df, ncp)
Rgeom.mat(nrow, ncol, prob)
Rpareto.mat(nrow, ncol, shape = 1, scale = 1)
Rfrechet.mat(nrow, ncol, lambda = 1, mu = 0, sigma = 1)
Rlaplace.mat(nrow, ncol, mu = 0, sigma = 1)
Rgumble.mat(nrow, ncol, mu = 0, sigma = 1)
Rarcsine.mat(nrow, ncol, min = 0, max = 1)
Rnorm.mat(nrow, ncol)
colRunif(nrow, ncol, min = rep_len(0, ncol), max = rep_len(1, ncol))
colRbeta(nrow, ncol, alpha, beta)
colRexp(nrow, ncol, rate = rep_len(1, ncol))
colRchisq(nrow, ncol, df)
colRgamma(nrow, ncol, shape, rate = rep_len(1, ncol))
colRgeom(nrow, ncol, prob)
colRcauchy(nrow, ncol, location = rep_len(0, ncol), scale = rep_len(1, ncol))
colRt(nrow, ncol, df, ncp)
colRpareto(nrow, ncol, shape = rep_len(1, ncol), scale = rep_len(1, ncol))
colRfrechet(nrow, ncol, lambda = rep_len(1, ncol), 
mu = rep_len(0, ncol), sigma = rep_len(1, ncol))
colRlaplace(nrow, ncol, mu = rep_len(0, ncol), sigma = rep_len(1, ncol))
colRgumble(nrow, ncol, mu = rep_len(0, ncol), sigma = rep_len(1, ncol))
colRarcsine(nrow, ncol, min = rep_len(0, ncol), max = rep_len(1, ncol))
colRnorm(nrow, ncol)
set_seed(seed = nanoTime())
}

\arguments{
\item{n}{
The number of values to generate.
}
\item{nrow}{
The number of rows.
}
\item{ncol}{
The number of columns.
}
\item{min}{
The lower value.
}
\item{max}{
The upper value.
}
\item{alpha}{
The shape parameter alpha.
}
\item{beta}{
The shape parameter beta.
}
\item{rate}{
\code{Rgamma}: The rate parameter.
}
\item{df}{
\code{Rt}: The degrees of freedom.
}
\item{lambda,shape}{
The shape parameter.
}
\item{location,mu}{
The location parameter.
}
\item{scale,sigma}{
The scale parameter.
}
\item{ncp}{
The non-centrality parameter.
}
\item{prob}{
The probability of success on each trial.
}
\item{seed}{
A single value, interpreted as an integer.
}
}

\details{
\itemize{
\item \code{Runif, Runif.mat, colRunif}: generates random values from the uniform distribution, similar to R's built-in \link{runif} function. The type used is \eqn{min + (max - min) \cdot U}, where \code{U} is a uniform random variable in the interval (0, 1).

\item \code{Rbeta, Rbeta.mat, colRbeta}: generates random values from the beta distribution with parameters \code{alpha} and \code{beta}. The type used involves generating two gamma-distributed variables \( X_1 \) and \( X_2 \), and returning \eqn{\frac{X_1}{X_1 + X_2}}, where \( X_1 \) and \( X_2 \) have the shape parameters \code{alpha} and \code{beta}, respectively.

\item \code{Rexp, Rexp.mat, colRexp}: generates random values from the exponential distribution with the specified \code{rate} parameter. The type used is \eqn{-\frac{\log(U)}{\text{rate}}}, where \code{U} is a uniform random variable in the interval (0, 1).

\item \code{Rchisq, Rchisq.mat, colRchisq}: generates random values from the chi-squared distribution with \code{df} degrees of freedom. The type used is the sum of the squares of \code{df} independent standard normal random variables, i.e., \eqn{Y_1^2 + Y_2^2 + \dots + Y_{df}^2}, where each \( Y_i \) is a standard normal random variable.

\item \code{Rgamma, Rgamma.mat, colRgamma}: generates random values from the gamma distribution with \code{shape} and \code{rate} parameters. The type used for shape greater than or equal to 1 is the Marsaglia and Tsang method, which involves generating a variable \( V \) and returning \eqn{d \cdot V \cdot \text{rate}}, where \( d \) is a function of the shape and \( V \) is derived from a normal random variable. For shape less than 1, a combination of the uniform and exponential distributions is used, involving an acceptance-rejection method.

\item \code{Rcauchy, Rcauchy.mat, colRcauchy}: generates random values from the Cauchy distribution with specified \code{location} and \code{scale} parameters. The type used for this is \eqn{\text{location} + \text{scale} \cdot \tan(\pi \cdot (U - 0.5))}, where \code{U} is a uniform random variable.

\item \code{Rt, Rt.mat, colRt}: generates random values from the t-distribution with specified \code{df} degrees of freedom and an optional non-centrality parameter \code{ncp}. The type used is \eqn{\frac{Z}{\sqrt{\frac{Y}{df}}}}, where \code{Z} is a standard normal random variable and \code{Y} is a chi-squared random variable. If \code{ncp} is provided, the type used is \eqn{\frac{Z + ncp}{\sqrt{\frac{Y}{df}}}}.

\item \code{Rgeom, Rgeom.mat, colRgeom}: generates random values from the geometric distribution with the specified probability \code{prob}. The type used is \eqn{\left\lfloor \frac{\log(U)}{\log(1 - prob)} \right\rfloor}, where \code{U} is a uniform random variable.

\item \code{Rpareto, Rpareto.mat, colRpareto}: generates random values from the Pareto distribution with parameters \code{shape} and \code{scale}. The type used is \eqn{\text{scale} \cdot (1 - U)^{-\frac{1}{\text{shape}}}}, where \code{U} is a uniform random variable in the interval (0, 1).

\item \code{Rfrechet, Rfrechet.mat, colRfrechet}: generates random values from the Frechet distribution with parameters \code{shape}, \code{mu}, and \code{scale}. The type used is \eqn{\mu + \text{scale} \cdot (-\log U)^{1/\text{shape}}}, where \code{U} is a uniform random variable in (0, 1).

\item \code{Rlaplace, Rlaplace.mat, colRlaplace}: generates random values from the Laplace distribution with location parameter \code{mu} and scale parameter \code{sigma}. The type used is \eqn{\mu - \sigma \cdot \text{sign}(U) \cdot \log(1 - 2 \cdot |U|)}, where \code{U} is a uniform random variable in (0, 1).

\item \code{Rgumble, Rgumble.mat, colRgumble}: generates random values from the Gumbel (type I extreme value) distribution with parameters \code{mu} and \code{sigma}. The type used is \eqn{\mu - \sigma \cdot \log(-\log U)}, where \code{U} is a uniform random variable in (0, 1).

\item \code{Rarcsine, Rarcsine.mat, colRarcsine}: generates random values from the arcsine distribution over the interval \code{[min, max]}. The type used is \eqn{\text{min} + (\text{max} - \text{min}) \cdot \sin^2(\pi U / 2)}, where \code{U} is a uniform random variable in (0, 1).

\item \code{set_seed}: Set the seed for the Rngs. Not working with parallel version of column - row sample.
}
}


\value{
Each \code{function}, for example \code{Runif}, returns a vector with simulated values from the respective distribution.
Each \code{function.mat}, for example \code{Runif.mat}, returns a matrix with simulated values from the respective distribution.
Each \code{colFunction}, for example \code{colRunif}, returns a matrix with column major simulated values from the respective distribution.
}

\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\link{runif}, \link{rbeta}, \link{rexp}, \link{rchisq}, \link{rgamma}, \link{rcauchy}, \link{rt}, \link{rgeom}
}

\examples{
# Scalar draws
x_unif     <- Runif(5, 0, 1)
x_beta     <- Rbeta(5, 2, 5)
x_exp      <- Rexp(5, 1.5)
x_chisq    <- Rchisq(5, 4)
x_gamma    <- Rgamma(5, 2, 2)
x_cauchy   <- Rcauchy(5, 0, 1)
x_t        <- Rt(5, df = 5, ncp = 2)
x_geom     <- Rgeom(5, 0.5)
x_pareto   <- Rpareto(5, shape = 2, scale = 1)
x_frechet  <- Rfrechet(5, lambda = 1, mu = 0, sigma = 1)
x_laplace  <- Rlaplace(5, mu = 0, sigma = 1)
x_gumblet  <- Rgumble(5, mu = 0, sigma = 1)
x_arcsine  <- Rarcsine(5, min = 0, max = 1)
x_norm     <- Rnorm(5)

#matrices
x_unif     <- Runif.mat(5,2, 0, 1)
x_beta     <- Rbeta.mat(5,2, 2, 5)
x_exp      <- Rexp.mat(5,2, 1.5)
x_chisq    <- Rchisq.mat(5,2, 4)
x_gamma    <- Rgamma.mat(5,2, 2, 2)
x_cauchy   <- Rcauchy.mat(5,2, 0, 1)
x_t        <- Rt.mat(5,2, df = 5, ncp = 2)
x_geom     <- Rgeom.mat(5,2, 0.5)
x_pareto   <- Rpareto.mat(5,2, shape = 2, scale = 1)
x_frechet  <- Rfrechet.mat(5,2, lambda = 1, mu = 0, sigma = 1)
x_laplace  <- Rlaplace.mat(5,2, mu = 0, sigma = 1)
x_gumblet  <- Rgumble.mat(5,2, mu = 0, sigma = 1)
x_arcsine  <- Rarcsine.mat(5,2, min = 0, max = 1)
x_norm     <- Rnorm.mat(5,2)

# Column-wise (vectorized by column) draws
x_col_unif    <- colRunif(5, 2, min = c(0, 1), max = c(1, 2))
x_col_beta    <- colRbeta(5, 2, alpha = c(2, 5), beta = c(5, 2))
x_col_exp     <- colRexp(5, 2, rate = c(1.5, 2.0))
x_col_chisq   <- colRchisq(5, 2, df = c(4, 5))
x_col_gamma   <- colRgamma(5, 2, shape = c(2, 3), rate = c(2, 1))
x_col_cauchy  <- colRcauchy(5, 2, location = c(0, 1), scale = c(1, 2))
x_col_t       <- colRt(5, 2, df = c(5, 6), ncp = c(2, 1))
x_col_geom    <- colRgeom(5, 2, prob = c(0.5, 0.3))
x_col_pareto  <- colRpareto(5, 2, shape = c(2, 3), scale = c(1, 1))
x_col_frechet <- colRfrechet(5, 2, lambda = c(1, 2), mu = c(0, 0), sigma = c(1, 1))
x_col_laplace <- colRlaplace(5, 2, mu = c(0, 1), sigma = c(1, 2))
x_col_gumble  <- colRgumble(5, 2, mu = c(0, 1), sigma = c(1, 2))
x_col_arcsine <- colRarcsine(5, 2, min = c(0, 0.5), max = c(1, 1.5))
x_col_norm    <- colRnorm(5, 2)

}
