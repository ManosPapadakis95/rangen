#include "ragen.h"
#include <Rcpp.h>

//[[Rcpp::export(name = "Runif")]]
Rcpp::NumericVector runif(size_t n, double min, double max) {
    return ragen::runif<Rcpp::NumericVector>(n, min, max);
}

//[[Rcpp::export(name = "Sample.int")]]
Rcpp::IntegerVector sample_int(size_t n, size_t size, bool replace) {
    return ragen::sample<Rcpp::IntegerVector>(n, size, replace);
}

//[[Rcpp::export(name = "Sample")]]
Rcpp::NumericVector sample(Rcpp::NumericVector x, size_t size, bool replace) {
    return ragen::sample<Rcpp::NumericVector>(x, size, replace);
}

//[[Rcpp::export(name = "Rbeta")]]
Rcpp::NumericVector rbeta(size_t size, double alpha, double beta) {
    return ragen::rbeta<Rcpp::NumericVector>(size, alpha, beta);
}

//[[Rcpp::export(name = "Rexp")]]
Rcpp::NumericVector rexp(size_t size, double rate) {
    return ragen::rexp<Rcpp::NumericVector>(size, rate);
}

//[[Rcpp::export(name = "Rchisq")]]
Rcpp::NumericVector rchisq(size_t size, double df) {
    return ragen::rchisq<Rcpp::NumericVector>(size, df);
}

//[[Rcpp::export(name = "Rgamma")]]
Rcpp::NumericVector rgamma(size_t size, double shape, double rate) {
    return ragen::rgamma<Rcpp::NumericVector>(size, shape, rate);
}

//[[Rcpp::export(name = "Rgeom")]]
Rcpp::NumericVector rgeom(size_t size, double prob) {
    return ragen::rgeom<Rcpp::NumericVector>(size, prob);
}

//[[Rcpp::export(name = "Rcauchy")]]
Rcpp::NumericVector rcauchy(size_t size, double location, double scale) {
    return ragen::rcauchy<Rcpp::NumericVector>(size, location, scale);
}

//[[Rcpp::export(name = "Rt")]]
Rcpp::NumericVector rt(size_t size, double df, double ncp) {
    return ragen::rt<Rcpp::NumericVector>(size, df, ncp);
}
