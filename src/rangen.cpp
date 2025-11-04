#include "rangen.h"
#include <RcppArmadillo.h>

//[[Rcpp::export(name = "Runif", signature = {n, min = 0, max = 1})]]
Rcpp::NumericVector runif(size_t n, double min = 0.0, double max = 1.0)
{
    return rangen::runif<Rcpp::NumericVector>(n, min, max);
}

//[[Rcpp::export(name = "Rbeta", signature = {n, alpha, beta})]]
Rcpp::NumericVector rbeta(size_t n, double alpha, double beta)
{
    return rangen::rbeta<Rcpp::NumericVector>(n, alpha, beta);
}

//[[Rcpp::export(name = "Rexp", signature = {n, rate = 1})]]
Rcpp::NumericVector rexp(size_t n, double rate = 1.0)
{
    return rangen::rexp<Rcpp::NumericVector>(n, rate);
}

//[[Rcpp::export(name = "Rchisq", signature = {n, df})]]
Rcpp::NumericVector rchisq(size_t n, double df)
{
    return rangen::rchisq<Rcpp::NumericVector>(n, df);
}

//[[Rcpp::export(name = "Rgamma", signature = {n, shape, rate = 1})]]
Rcpp::NumericVector rgamma(size_t n, double shape, double rate = 1.0)
{
    return rangen::rgamma<Rcpp::NumericVector>(n, shape, rate);
}

//[[Rcpp::export(name = "Rgeom", signature = {n, prob})]]
Rcpp::NumericVector rgeom(size_t n, double prob)
{
    return rangen::rgeom<Rcpp::NumericVector>(n, prob);
}

//[[Rcpp::export(name = "Rcauchy", signature = {n, location = 0, scale = 1})]]
Rcpp::NumericVector rcauchy(size_t n, double location = 0.0, double scale = 1.0)
{
    return rangen::rcauchy<Rcpp::NumericVector>(n, location, scale);
}

//[[Rcpp::export(name = "Rt", signature = {n, df, ncp})]]
Rcpp::NumericVector rt(size_t n, double df, double ncp)
{
    return rangen::rt<Rcpp::NumericVector>(n, df, ncp);
}

//[[Rcpp::export(name = "Rpareto", signature = {n, shape = 1, scale = 1})]]
Rcpp::NumericVector rpareto(size_t n, double shape = 1.0, double scale = 1.0)
{
    return rangen::rpareto<Rcpp::NumericVector>(n, shape, scale);
}

//[[Rcpp::export(name = "Rfrechet", signature = {n, lambda = 1, mu = 0, sigma = 1})]]
Rcpp::NumericVector rfrechet(size_t n, double lambda = 1.0, double mu = 0.0, double sigma = 1.0)
{
    return rangen::rfrechet<Rcpp::NumericVector>(n, lambda, mu, sigma);
}

//[[Rcpp::export(name = "Rlaplace", signature = {n, mu = 0, sigma = 1})]]
Rcpp::NumericVector rlaplace(size_t n, double mu = 0.0, double sigma = 1.0)
{
    return rangen::rlaplace<Rcpp::NumericVector>(n, mu, sigma);
}

//[[Rcpp::export(name = "Rgumble", signature = {n, mu = 0, sigma = 1})]]
Rcpp::NumericVector rgumble(size_t n, double mu = 0.0, double sigma = 1.0)
{
    return rangen::rgumble<Rcpp::NumericVector>(n, mu, sigma);
}

//[[Rcpp::export(name = "Rarcsine", signature = {n, min = 0, max = 1})]]
Rcpp::NumericVector rarcsine(size_t n, double min = 0.0, double max = 1.0)
{
    return rangen::rarcsine<Rcpp::NumericVector>(n, min, max);
}

//[[Rcpp::export(name = "Rnorm")]]
Rcpp::NumericVector rnorm(size_t n)
{
    return rangen::rnorm<Rcpp::NumericVector>(n);
}

//[[Rcpp::export(signature = {nrow, ncol, min = rep_len(0, ncol), max = rep_len(1, ncol)})]]
Rcpp::NumericVector colRunif(size_t nrow, size_t ncol, Rcpp::NumericVector min, Rcpp::NumericVector max)
{
    return rangen::colRng(nrow, ncol, rangen::runif<arma::colvec>, min, max);
}

//[[Rcpp::export(signature = {nrow, ncol, alpha, beta})]]
Rcpp::NumericVector colRbeta(size_t nrow, size_t ncol, Rcpp::NumericVector alpha, Rcpp::NumericVector beta)
{
    return rangen::colRng(nrow, ncol, rangen::rbeta<arma::colvec>, alpha, beta);
}

//[[Rcpp::export(signature = {nrow, ncol, rate = rep_len(1, ncol)})]]
Rcpp::NumericVector colRexp(size_t nrow, size_t ncol, Rcpp::NumericVector rate)
{
    return rangen::colRng(nrow, ncol, rangen::rexp<arma::colvec>, rate);
}

//[[Rcpp::export(signature = {nrow, ncol, df})]]
Rcpp::NumericVector colRchisq(size_t nrow, size_t ncol, Rcpp::NumericVector df)
{
    return rangen::colRng(nrow, ncol, rangen::rchisq<arma::colvec>, df);
}

//[[Rcpp::export(signature = {nrow, ncol, shape, rate = rep_len(1, ncol)})]]
Rcpp::NumericVector colRgamma(size_t nrow, size_t ncol, Rcpp::NumericVector shape, Rcpp::NumericVector rate)
{
    return rangen::colRng(nrow, ncol, rangen::rgamma<arma::colvec>, shape, rate);
}

//[[Rcpp::export(signature = {nrow, ncol, prob})]]
Rcpp::NumericVector colRgeom(size_t nrow, size_t ncol, Rcpp::NumericVector prob)
{
    return rangen::colRng(nrow, ncol, rangen::rgeom<arma::colvec>, prob);
}

//[[Rcpp::export(signature = {nrow, ncol, location = rep_len(0, ncol), scale = rep_len(1, ncol)})]]
Rcpp::NumericVector colRcauchy(size_t nrow, size_t ncol, Rcpp::NumericVector location, Rcpp::NumericVector scale)
{
    return rangen::colRng(nrow, ncol, rangen::rcauchy<arma::colvec>, location, scale);
}

//[[Rcpp::export(signature = {nrow, ncol, df, ncp})]]
Rcpp::NumericVector colRt(size_t nrow, size_t ncol, Rcpp::NumericVector df, Rcpp::NumericVector ncp)
{
    return rangen::colRng(nrow, ncol, rangen::rt<arma::colvec>, df, ncp);
}

//[[Rcpp::export(signature = {nrow, ncol, shape = rep_len(1, ncol), scale = rep_len(1, ncol)})]]
Rcpp::NumericVector colRpareto(size_t nrow, size_t ncol, Rcpp::NumericVector shape, Rcpp::NumericVector scale)
{
    return rangen::colRng(nrow, ncol, rangen::rpareto<arma::colvec>, shape, scale);
}

//[[Rcpp::export(signature = {nrow, ncol, lambda = rep_len(1, ncol), mu = rep_len(0, ncol), sigma = rep_len(1, ncol)})]]
Rcpp::NumericVector colRfrechet(size_t nrow, size_t ncol, Rcpp::NumericVector lambda, Rcpp::NumericVector mu, Rcpp::NumericVector sigma)
{
    return rangen::colRng(nrow, ncol, rangen::rfrechet<arma::colvec>, lambda, mu, sigma);
}

//[[Rcpp::export(signature = {nrow, ncol, mu = rep_len(0, ncol), sigma = rep_len(1, ncol)})]]
Rcpp::NumericVector colRlaplace(size_t nrow, size_t ncol, Rcpp::NumericVector mu, Rcpp::NumericVector sigma)
{
    return rangen::colRng(nrow, ncol, rangen::rlaplace<arma::colvec>, mu, sigma);
}

//[[Rcpp::export(signature = {nrow, ncol, mu = rep_len(0, ncol), sigma = rep_len(1, ncol)})]]
Rcpp::NumericVector colRgumble(size_t nrow, size_t ncol, Rcpp::NumericVector mu, Rcpp::NumericVector sigma)
{
    return rangen::colRng(nrow, ncol, rangen::rgumble<arma::colvec>, mu, sigma);
}

//[[Rcpp::export(signature = {nrow, ncol, min = rep_len(0, ncol), max = rep_len(1, ncol)})]]
Rcpp::NumericVector colRarcsine(size_t nrow, size_t ncol, Rcpp::NumericVector min, Rcpp::NumericVector max)
{
    return rangen::colRng(nrow, ncol, rangen::rarcsine<arma::colvec>, min, max);
}

//[[Rcpp::export]]
Rcpp::NumericVector colRnorm(size_t nrow, size_t ncol)
{
    return rangen::colRng(nrow, ncol, rangen::rnorm<arma::colvec>);
}

//[[Rcpp::export(name = "Runif.mat", signature = {nrow, ncol, min = 0, max = 1})]]
Rcpp::NumericMatrix runif_mat(size_t nrow, size_t ncol, double min = 0.0, double max = 1.0)
{
    return rangen::runif_mat<Rcpp::NumericMatrix>(nrow, ncol, min, max);
}

//[[Rcpp::export(name = "Rbeta.mat", signature = {nrow, ncol, alpha, beta})]]
Rcpp::NumericMatrix rbeta_mat(size_t nrow, size_t ncol, double alpha, double beta) {
    return rangen::rbeta_mat<Rcpp::NumericMatrix>(nrow, ncol, alpha, beta);
}

//[[Rcpp::export(name = "Rexp.mat", signature = {nrow, ncol, rate = 1})]]
Rcpp::NumericMatrix rexp_mat(size_t nrow, size_t ncol, double rate = 1.0) {
    return rangen::rexp_mat<Rcpp::NumericMatrix>(nrow, ncol, rate);
}

//[[Rcpp::export(name = "Rchisq.mat", signature = {nrow, ncol, df})]]
Rcpp::NumericMatrix rchisq_mat(size_t nrow, size_t ncol, double df) {
    return rangen::rchisq_mat<Rcpp::NumericMatrix>(nrow, ncol, df);
}

//[[Rcpp::export(name = "Rgamma.mat", signature = {nrow, ncol, shape, rate = 1})]]
Rcpp::NumericMatrix rgamma_mat(size_t nrow, size_t ncol, double shape, double rate = 1.0) {
    return rangen::rgamma_mat<Rcpp::NumericMatrix>(nrow, ncol, shape, rate);
}

//[[Rcpp::export(name = "Rgeom.mat", signature = {nrow, ncol, prob})]]
Rcpp::NumericMatrix rgeom_mat(size_t nrow, size_t ncol, double prob) {
    return rangen::rgeom_mat<Rcpp::NumericMatrix>(nrow, ncol, prob);
}

//[[Rcpp::export(name = "Rcauchy.mat", signature = {nrow, ncol, location = 0, scale = 1})]]
Rcpp::NumericMatrix rcauchy_mat(size_t nrow, size_t ncol, double location = 0.0, double scale = 1.0) {
    return rangen::rcauchy_mat<Rcpp::NumericMatrix>(nrow, ncol, location, scale);
}

//[[Rcpp::export(name = "Rt.mat", signature = {nrow, ncol, df, ncp})]]
Rcpp::NumericMatrix rt_mat(size_t nrow, size_t ncol, double df, double ncp) {
    return rangen::rt_mat<Rcpp::NumericMatrix>(nrow, ncol, df, ncp);
}

//[[Rcpp::export(name = "Rpareto.mat", signature = {nrow, ncol, shape = 1, scale = 1})]]
Rcpp::NumericMatrix rpareto_mat(size_t nrow, size_t ncol, double shape = 1.0, double scale = 1.0) {
    return rangen::rpareto_mat<Rcpp::NumericMatrix>(nrow, ncol, shape, scale);
}

//[[Rcpp::export(name = "Rfrechet.mat", signature = {nrow, ncol, lambda = 1, mu = 0, sigma = 1})]]
Rcpp::NumericMatrix rfrechet_mat(size_t nrow, size_t ncol, double lambda = 1.0, double mu = 0.0, double sigma = 1.0) {
    return rangen::rfrechet_mat<Rcpp::NumericMatrix>(nrow, ncol, lambda, mu, sigma);
}

//[[Rcpp::export(name = "Rlaplace.mat", signature = {nrow, ncol, mu = 0, sigma = 1})]]
Rcpp::NumericMatrix rlaplace_mat(size_t nrow, size_t ncol, double mu = 0.0, double sigma = 1.0) {
    return rangen::rlaplace_mat<Rcpp::NumericMatrix>(nrow, ncol, mu, sigma);
}

//[[Rcpp::export(name = "Rgumble.mat", signature = {nrow, ncol, mu = 0, sigma = 1})]]
Rcpp::NumericMatrix rgumble_mat(size_t nrow, size_t ncol, double mu = 0.0, double sigma = 1.0) {
    return rangen::rgumble_mat<Rcpp::NumericMatrix>(nrow, ncol, mu, sigma);
}

//[[Rcpp::export(name = "Rarcsine.mat", signature = {nrow, ncol, min = 0, max = 1})]]
Rcpp::NumericMatrix rarcsine_mat(size_t nrow, size_t ncol, double min = 0.0, double max = 1.0) {
    return rangen::rarcsine_mat<Rcpp::NumericMatrix>(nrow, ncol, min, max);
}

//[[Rcpp::export(name = "Rnorm.mat")]]
Rcpp::NumericMatrix rnorm_mat(size_t nrow, size_t ncol) {
    return rangen::rnorm_mat<Rcpp::NumericMatrix>(nrow, ncol);
}

//[[Rcpp::export(name = "Sample.int", signature = {n, size = n, replace = FALSE})]]
Rcpp::IntegerVector sample_int(size_t n, size_t size, const bool replace = false)
{
    return rangen::sample<Rcpp::IntegerVector>(n, size, replace);
}

//[[Rcpp::export(name = "Sample", signature = {x, size = length(x), replace = FALSE})]]
Rcpp::NumericVector sample(Rcpp::NumericVector x, size_t size, const bool replace = false)
{
    return rangen::sample<Rcpp::NumericVector>(x, size, replace);
}

//[[Rcpp::export(name = "colSample", signature = {x, size = rep_len(nrow(x), ncol(x)), replace = rep_len(FALSE, ncol(x)), parallel = FALSE, cores = 0})]]
Rcpp::NumericMatrix colSample(Rcpp::NumericMatrix x, Rcpp::NumericVector size, Rcpp::LogicalVector replace, const bool parallel = false, const size_t cores = 0)
{
    return rangen::colSample(x, size, replace, parallel, cores);
}

//[[Rcpp::export(name = "rowSample", signature = {x, size = rep_len(ncol(x), nrow(x)), replace = rep_len(FALSE, nrow(x)), parallel = FALSE, cores = 0})]]
Rcpp::NumericMatrix rowSample(Rcpp::NumericMatrix x, Rcpp::NumericVector size, Rcpp::LogicalVector replace, const bool parallel = false, const size_t cores = 0)
{
    return rangen::rowSample(x, size, replace, parallel, cores);
}

//[[Rcpp::export]]
double nanoTime()
{
    return rangen::internal::get_cur_nano();
}

//[[Rcpp::export(signature = {seed = nanoTime()})]]
void set_seed(double seed)
{
    return rangen::setSeed(static_cast<size_t>(seed));
}