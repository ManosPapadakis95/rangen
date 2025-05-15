#include "rangen.h"
#include <Rcpp.h>

//[[Rcpp::export(name = "Runif", signature = {n, min = 0, max = 1})]]
Rcpp::NumericVector runif(size_t n, double min = 0.0, double max = 1.0) {
    return rangen::runif<Rcpp::NumericVector>(n, min, max);
}

//[[Rcpp::export(name = "Rbeta", signature = {n, alpha, beta})]]
Rcpp::NumericVector rbeta(size_t n, double alpha, double beta) {
    return rangen::rbeta<Rcpp::NumericVector>(n, alpha, beta);
}

//[[Rcpp::export(name = "Rexp", signature = {n, rate = 1})]]
Rcpp::NumericVector rexp(size_t n, double rate = 1.0) {
    return rangen::rexp<Rcpp::NumericVector>(n, rate);
}

//[[Rcpp::export(name = "Rchisq", signature = {n, df})]]
Rcpp::NumericVector rchisq(size_t n, double df) {
    return rangen::rchisq<Rcpp::NumericVector>(n, df);
}

//[[Rcpp::export(name = "Rgamma", signature = {n, shape, rate = 1})]]
Rcpp::NumericVector rgamma(size_t n, double shape, double rate = 1.0) {
    return rangen::rgamma<Rcpp::NumericVector>(n, shape, rate);
}

//[[Rcpp::export(name = "Rgeom", signature = {n, prob})]]
Rcpp::NumericVector rgeom(size_t n, double prob) {
    return rangen::rgeom<Rcpp::NumericVector>(n, prob);
}

//[[Rcpp::export(name = "Rcauchy", signature = {n, location = 0, scale = 1})]]
Rcpp::NumericVector rcauchy(size_t n, double location = 0.0, double scale = 1.0) {
    return rangen::rcauchy<Rcpp::NumericVector>(n, location, scale);
}

//[[Rcpp::export(name = "Rt", signature = {n, df, ncp})]]
Rcpp::NumericVector rt(size_t n, double df, double ncp) {
    return rangen::rt<Rcpp::NumericVector>(n, df, ncp);
}

//[[Rcpp::export(name = "Rpareto", signature = {n, shape = 1, scale = 1})]]
Rcpp::NumericVector rpareto(size_t n, double shape = 1.0, double scale = 1.0) {
    return rangen::rpareto<Rcpp::NumericVector>(n, shape, scale);
}

//[[Rcpp::export(name = "Rfrechet", signature = {n, lambda = 1, mu = 0, sigma = 1})]]
Rcpp::NumericVector rfrechet(size_t n, double lambda = 1.0, double mu = 0.0, double sigma = 1.0) {
    return rangen::rfrechet<Rcpp::NumericVector>(n, lambda, mu, sigma);
}

//[[Rcpp::export(name = "Rlaplace", signature = {n, mu = 0, sigma = 1})]]
Rcpp::NumericVector rlaplace(size_t n, double mu = 0.0, double sigma = 1.0) {
    return rangen::rlaplace<Rcpp::NumericVector>(n, mu, sigma);
}

//[[Rcpp::export(name = "Rgumblet", signature = {n, mu = 0, sigma = 1})]]
Rcpp::NumericVector rgumble(size_t n, double mu = 0.0, double sigma = 1.0) {
    return rangen::rgumble<Rcpp::NumericVector>(n, mu, sigma);
}

//[[Rcpp::export(name = "Rarcsine", signature = {n, min = 0, max = 1})]]
Rcpp::NumericVector rarcsine(size_t n, double min = 0.0, double max = 1.0) {
    return rangen::rarcsine<Rcpp::NumericVector>(n, min, max);
}

//[[Rcpp::export(name = "Sample.int", signature = {n, size = n, replace = FALSE})]]
Rcpp::IntegerVector sample_int(size_t n, size_t size, bool replace = false) {
    return rangen::sample<Rcpp::IntegerVector>(n, size, replace);
}

//[[Rcpp::export(name = "Sample", signature = {x, size = n, replace = FALSE})]]
Rcpp::NumericVector sample(Rcpp::NumericVector x, size_t size, bool replace = false) {
    return rangen::sample<Rcpp::NumericVector>(x, size, replace);
}

//[[Rcpp::export(name = "colSample", signature = {x, size = rep_len(n, ncol(x)), replace = rep_len(FALSE, ncol(x))})]]
Rcpp::NumericMatrix colSample(Rcpp::NumericMatrix x, Rcpp::NumericMatrix size, Rcpp::LogicalVector replace) {
    return rangen::colSample(x,size,replace);
}

//[[Rcpp::export(name = "rowSample", signature = {x, size = rep_len(n, nrow(x)), replace = rep_len(FALSE, nrow(x))})]]
Rcpp::NumericMatrix rowSample(Rcpp::NumericMatrix x, Rcpp::NumericMatrix size, Rcpp::LogicalVector replace) {
    return rangen::rowSample(x,size,replace);
}