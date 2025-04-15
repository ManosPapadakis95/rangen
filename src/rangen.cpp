#include "rangen.h"
#include <Rcpp.h>

//[[Rcpp::export(name = "Runif")]]
Rcpp::NumericVector runif(size_t n, double min = 0.0, double max = 1.0) {
    return rangen::runif<Rcpp::NumericVector>(n, min, max);
}

//[[Rcpp::export(name = "Rbeta")]]
Rcpp::NumericVector rbeta(size_t size, double alpha, double beta) {
    return rangen::rbeta<Rcpp::NumericVector>(size, alpha, beta);
}

//[[Rcpp::export(name = "Rexp")]]
Rcpp::NumericVector rexp(size_t size, double rate = 1.0) {
    return rangen::rexp<Rcpp::NumericVector>(size, rate);
}

//[[Rcpp::export(name = "Rchisq")]]
Rcpp::NumericVector rchisq(size_t size, double df) {
    return rangen::rchisq<Rcpp::NumericVector>(size, df);
}

//[[Rcpp::export(name = "Rgamma")]]
Rcpp::NumericVector rgamma(size_t size, double shape, double rate = 1.0) {
    return rangen::rgamma<Rcpp::NumericVector>(size, shape, rate);
}

//[[Rcpp::export(name = "Rgeom")]]
Rcpp::NumericVector rgeom(size_t size, double prob) {
    return rangen::rgeom<Rcpp::NumericVector>(size, prob);
}

//[[Rcpp::export(name = "Rcauchy")]]
Rcpp::NumericVector rcauchy(size_t size, double location = 0.0, double scale = 1.0) {
    return rangen::rcauchy<Rcpp::NumericVector>(size, location, scale);
}

//[[Rcpp::export(name = "Rt")]]
Rcpp::NumericVector rt(size_t size, double df, double ncp) {
    return rangen::rt<Rcpp::NumericVector>(size, df, ncp);
}

//[[Rcpp::export(name = "Rpareto")]]
Rcpp::NumericVector rpareto(size_t size, double shape = 1.0, double scale = 1.0) {
    return rangen::rpareto<Rcpp::NumericVector>(size, shape, scale);
}

//[[Rcpp::export(name = "Rfrechet")]]
Rcpp::NumericVector rfrechet(size_t size, double shape = 1.0, double mean = 0.0, double scale = 1.0) {
    return rangen::rfrechet<Rcpp::NumericVector>(size, shape, mean, scale);
}

//[[Rcpp::export(name = "Rlaplace")]]
Rcpp::NumericVector rlaplace(size_t size, double mean = 0.0, double sigma = 1.0) {
    return rangen::rlaplace<Rcpp::NumericVector>(size, mean, sigma);
}

//[[Rcpp::export(name = "Rgumblet")]]
Rcpp::NumericVector rgumble(size_t size, double mean = 0.0, double sigma = 1.0) {
    return rangen::rgumble<Rcpp::NumericVector>(size, mean, sigma);
}

//[[Rcpp::export(name = "Rarcsine")]]
Rcpp::NumericVector rarcsine(size_t size, double min = 0.0, double max = 1.0) {
    return rangen::rarcsine<Rcpp::NumericVector>(size, min, max);
}

//[[Rcpp::export(name = "Sample.int")]]
Rcpp::IntegerVector sample_int(size_t n, size_t size, bool replace = false) {
    return rangen::sample<Rcpp::IntegerVector>(n, size, replace);
}

//[[Rcpp::export(name = "Sample")]]
Rcpp::NumericVector sample(Rcpp::NumericVector x, size_t size, bool replace = false) {
    return rangen::sample<Rcpp::NumericVector>(x, size, replace);
}

//[[Rcpp::export(name = "colSample", signature = {x, size, replace = rep_len(FALSE, ncol(x))})]]
Rcpp::NumericVector colSample(Rcpp::NumericMatrix x, size_t size, Rcpp::LogicalVector replace) {
    return rangen::colSample(x,size,replace);
}

//[[Rcpp::export(name = "rowSample", signature = {x, size, replace = rep_len(FALSE, nrow(x))})]]
Rcpp::NumericVector rowSample(Rcpp::NumericMatrix x, size_t size, Rcpp::LogicalVector replace) {
    return rangen::rowSample(x,size,replace);
}