#pragma once

#include "distributions.h"
#include "assertions.hpp"

namespace rangen
{

    inline Rcpp::NumericMatrix colSample(Rcpp::NumericMatrix x, Rcpp::NumericVector size, Rcpp::LogicalVector replace)
    {
        const size_t n = rangen_internal::ncol(x);
        arma::mat xx = rangen_internal::getArmaFrom(x, false);
        Rcpp::NumericMatrix res(xx.n_rows, xx.n_cols);
        arma::mat Res = rangen_internal::getArmaFrom(res, false);

        for (size_t i = 0; i < n; ++i)
        {
            Res.col(i) = sample<arma::colvec>(xx.col(i), (size_t)size[i], replace[i]);
        }
        return res;
    }

    // Rcpp::List colSampleL(Rcpp::NumericMatrix x, Rcpp::NumericVector size, Rcpp::LogicalVector replace)
    // {
    //     const size_t n = rangen_internal::ncol(x);
    //     arma::mat xx = rangen_internal::getArmaFrom(x, false);
    //     Rcpp::List res(n);

    //     for (size_t i = 0; i < n; ++i)
    //     {
    //         Res[i] = sample<arma::colvec>(xx.col(i), (size_t)size[i], replace[i]);
    //     }
    //     return res;
    // }

    inline Rcpp::NumericMatrix rowSample(Rcpp::NumericMatrix x, Rcpp::NumericVector size, Rcpp::LogicalVector replace)
    {
        const size_t m = rangen_internal::nrow(x);
        arma::mat xx = rangen_internal::getArmaFrom(x, false);
        Rcpp::NumericMatrix res(xx.n_rows, xx.n_cols);
        arma::mat Res = rangen_internal::getArmaFrom(res, false);

        for (size_t i = 0; i < m; ++i)
        {
            Res.row(i) = sample<arma::rowvec>(xx.row(i), (size_t)size[i], replace[i]);
        }
        return res;
    }

    // Rcpp::List rowSample(Rcpp::NumericMatrix x, Rcpp::NumericVector size, Rcpp::LogicalVector replace)
    // {
    //     const size_t m = rangen_internal::nrow(x);
    //     arma::mat xx = rangen_internal::getArmaFrom(x, false);
    //     Rcpp::List res(m);

    //     for (size_t i = 0; i < m; ++i)
    //     {
    //         Res[i] = sample<arma::rowvec>(xx.row(i), (size_t)size[i], replace[i]);
    //     }
    //     return res;
    // }

}