\name{Random integer values simulation}
\alias{Sample.int}
\alias{Sample}
\alias{colSample}
\alias{rowSample}
\title{
Random integer values simulation
}

\description{
Random integer values simulation.
}

\usage{
Sample.int(n, size = n, replace = FALSE) 
Sample(x, size = n, replace = FALSE)
colSample(x, size = rep_len(n, ncol(x)), replace = rep_len(FALSE, ncol(x)), parallel = FALSE, cores = 0)
rowSample(x, size = rep_len(n, nrow(x)), replace = rep_len(FALSE, nrow(x)), parallel = FALSE, cores = 0)
}

\arguments{
\item{x}{
A numeric vector for sampling or a matrix for column-row sampling.
}
\item{n}{
This must be an integer value. The function will then draw random integer values from 1:n. 
}
\item{size}{
The number of integer values to sample.
}
\item{replace}{
Do you want to sample with replacement? If yes, set this equal to TRUE.
}
\item{parallel}{
Do you want to do it in parallel, for column - row major, in C++? TRUE or FALSE.
}
\item{cores}{
Number of cores to use for parallelism. Valid only when argument parallel is set to TRUE. 
Default value is 0 and it means the maximum supported cores.
}
}

\details{
These functions provide flexible sampling utilities similar in purpose to R's base functions \link{sample.int} and \link{sample}. Each function operates on a different structure:

\itemize{
  \item \code{Sample.int}: Generates a random sample of integers from \code{1} to \code{n}.
  \item \code{Sample}:     Samples elements from the vector \code{x}.
  \item \code{colSample}:  Performs column-wise sampling on a matrix or data frame, selecting \code{size[i]} elements from each column \code{i}.
  \item \code{rowSample}:  Performs row-wise sampling on a matrix or data frame, selecting \code{size[i]} elements from each row \code{i}.
}

All functions support sampling with or without replacement.
Parallel versions do not support seeding.
}

\value{
\item{Sample.int, Sample}{
A vector of sampled values.
}
\item{colSample, rowSample}{
A matrix or data frame containing the sampled elements.
}
}


\author{
R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{sample}, \link{sample.int} 
}
}

\examples{
# Sample integers from 1 to 10 with replacement (faster than base::sample.int)
x <- Sample.int(10, 1000, replace = TRUE)

# Sample from the vector 'x' (faster than base::sample)
xs <- Sample(x)

# Create a matrix and perform column-wise sampling
# 'size' must have the same length as the number of columns
mat <- matrix(1:20, nrow = 5, ncol = 4)
colSample(mat, size = rep(5, ncol(mat)), replace = TRUE)

# Create a matrix and perform row-wise sampling
# 'size' must have the same length as the number of rows
rowSample(mat, size = rep(4, nrow(mat)), replace = FALSE)
}

