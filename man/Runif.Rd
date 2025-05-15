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
\alias{Rgumblet}
\alias{Rarcsine}
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
Rgumblet(n, mu = 0, sigma = 1)
Rarcsine(n, min = 0, max = 1)
}

\arguments{
\item{n}{
The number of values to generate.
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
}

\details{
    
\itemize{
\item \code{Runif}: generates random values from the uniform distribution, similar to R's built-in \code{\link{runif}} function. The type used is \eqn{min + (max - min) \cdot U}, where \code{U} is a uniform random variable in the interval (0, 1).

\item \code{Rbeta}: generates random values from the beta distribution with parameters \code{alpha} and \code{beta}. The type used involves generating two gamma-distributed variables \( X_1 \) and \( X_2 \), and returning \eqn{\frac{X_1}{X_1 + X_2}}, where \( X_1 \) and \( X_2 \) have the shape parameters \code{alpha} and \code{beta}, respectively.

\item \code{Rexp}: generates random values from the exponential distribution with the specified \code{rate} parameter. The type used is \eqn{-\frac{\log(U)}{\text{rate}}}, where \code{U} is a uniform random variable in the interval (0, 1).

\item \code{Rchisq}: generates random values from the chi-squared distribution with \code{df} degrees of freedom. The type used is the sum of the squares of \code{df} independent standard normal random variables, i.e., \eqn{Y_1^2 + Y_2^2 + \dots + Y_{df}^2}, where each \( Y_i \) is a standard normal random variable.

\item \code{Rgamma}: generates random values from the gamma distribution with \code{shape} and \code{rate} parameters. The type used for shape greater than or equal to 1 is the Marsaglia and Tsang method, which involves generating a variable \( V \) and returning \eqn{d \cdot V \cdot \text{rate}}, where \( d \) is a function of the shape and \( V \) is derived from a normal random variable. For shape less than 1, a combination of the uniform and exponential distributions is used, involving an acceptance-rejection method.

\item \code{Rcauchy}: generates random values from the Cauchy distribution with specified \code{location} and \code{scale} parameters. The type used for this is \eqn{\text{location} + \text{scale} \cdot \tan(\pi \cdot (U - 0.5))}, where \code{U} is a uniform random variable.

\item \code{Rt}: generates random values from the t-distribution with specified \code{df} degrees of freedom and an optional non-centrality parameter \code{ncp}. The type used is \eqn{\frac{Z}{\sqrt{\frac{Y}{df}}}}, where \code{Z} is a standard normal random variable and \code{Y} is a chi-squared random variable. If \code{ncp} is provided, the type used is \eqn{\frac{Z + ncp}{\sqrt{\frac{Y}{df}}}}.

\item \code{Rgeom}: generates random values from the geometric distribution with the specified probability \code{prob}. The type used is \eqn{\left\lfloor \frac{\log(U)}{\log(1 - prob)} \right\rfloor}, where \code{U} is a uniform random variable.

\item \code{Rpareto}: generates random values from the Pareto distribution with parameters \code{shape} and \code{scale}. The type used is \eqn{\text{scale} \cdot (1 - U)^{-\frac{1}{\text{shape}}}}, where \code{U} is a uniform random variable in the interval (0, 1).

\item \code{Rfrechet}: generates random values from the Frechet distribution with parameters \code{shape}, \code{mu}, and \code{scale}. The type used is \eqn{\mu + \text{scale} \cdot (-\log U)^{1/\text{shape}}}, where \code{U} is a uniform random variable in (0, 1).

\item \code{Rlaplace}: generates random values from the Laplace distribution with location parameter \code{mu} and scale parameter \code{sigma}. The type used is \eqn{\mu - \sigma \cdot \text{sign}(U) \cdot \log(1 - 2 \cdot |U|)}, where \code{U} is a uniform random variable in (0, 1).

\item \code{Rgumblet}: generates random values from the Gumbel (type I extreme value) distribution with parameters \code{mu} and \code{sigma}. The type used is \eqn{\mu - \sigma \cdot \log(-\log U)}, where \code{U} is a uniform random variable in (0, 1).

\item \code{Rarcsine}: generates random values from the arcsine distribution over the interval \code{[min, max]}. The type used is \eqn{\text{min} + (\text{max} - \text{min}) \cdot \sin^2(\pi U / 2)}, where \code{U} is a uniform random variable in (0, 1).
}

}


\value{
Each function returns a vector with simulated values from the respective distribution.
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{runif}}, \code{\link{rbeta}}, \code{\link{rexp}}, \code{\link{rchisq}}, \code{\link{rgamma}}, \code{\link{rcauchy}}, \code{\link{rt}}, \code{\link{rgeom}}
}

\examples{
x_unif <- Runif(1000, 0, 1)
x_beta <- Rbeta(1000, 2, 5)
x_exp <- Rexp(1000, 1.5)
x_chisq <- Rchisq(1000, 4)
x_gamma <- Rgamma(1000, 2, 2)
x_cauchy <- Rcauchy(1000, 0, 1)
x_t <- Rt(1000, 5, 2)
x_geom <- Rgeom(1000, 0.5)
}
