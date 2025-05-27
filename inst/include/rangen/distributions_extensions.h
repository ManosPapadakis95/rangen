#pragma once

#include <RcppArmadillo.h>
#include "distributions.h"
#include "assertions.hpp"

namespace rangen
{

    namespace rangen_internal
    {
        template <class T>
        arma::Mat<typename std::remove_reference<typename T::value_type>::type> getArmaFrom(T x, const bool copy = true)
        {

            Assertion::has_subscript_operator<T>::check_concept();
            Assertion::has_size<T>::check_concept();

            arma::Mat<typename std::remove_reference<typename T::value_type>::type> res(x.begin(), nrow(x), ncol(x), copy);
            return res;
        }

        template <class T>
        arma::Col<typename std::remove_reference<typename T::value_type>::type> getArmaColFrom(T x, const bool copy = true)
        {

            Assertion::has_subscript_operator<T>::check_concept();
            Assertion::has_size<T>::check_concept();

            arma::Col<typename std::remove_reference<typename T::value_type>::type> res(x.begin(), size(x), copy);
            return res;
        }
    }

    inline Rcpp::NumericMatrix colSample(Rcpp::NumericMatrix x, Rcpp::NumericVector size, Rcpp::LogicalVector replace, const bool parallel = false, const size_t cores = rangen_internal::get_num_of_threads())
    {
        const size_t n = rangen_internal::ncol(x);
        const size_t m = *std::max_element(size.begin(), size.end());
        arma::mat xx = rangen_internal::getArmaFrom(x, false);
        arma::colvec ss = rangen_internal::getArmaColFrom(size, false);
        Rcpp::NumericMatrix res(m, n);
        arma::mat Res = rangen_internal::getArmaFrom(res, false);

        if (parallel)
        {

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
            for (size_t i = 0; i < n; ++i)
            {
                Res.col(i) = sample<arma::colvec, decltype(xx.col(i))>(xx.col(i), (size_t)ss[i], replace[i], parallel);
            }
        }
        else
        {
            for (size_t i = 0; i < n; ++i)
            {
                Res.col(i) = sample<arma::colvec, decltype(xx.col(i))>(xx.col(i), (size_t)ss[i], replace[i], parallel);
            }
        }
        return res;
    }

    inline Rcpp::NumericMatrix rowSample(Rcpp::NumericMatrix x, Rcpp::NumericVector size, Rcpp::LogicalVector replace, const bool parallel = false, const size_t cores = rangen_internal::get_num_of_threads())
    {
        const size_t m = rangen_internal::nrow(x);
        const size_t n = *std::max_element(size.begin(), size.end());
        arma::mat xx = rangen_internal::getArmaFrom(x, false);
        arma::colvec ss = rangen_internal::getArmaColFrom(size, false);
        Rcpp::NumericMatrix res(m, n);
        arma::mat Res = rangen_internal::getArmaFrom(res, false);

        if (parallel)
        {

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
            for (size_t i = 0; i < m; ++i)
            {
                Res.row(i) = sample<arma::rowvec, decltype(xx.row(i))>(xx.row(i), (size_t)ss[i], replace[i], parallel);
            }
        }
        else
        {
            for (size_t i = 0; i < m; ++i)
            {
                Res.row(i) = sample<arma::rowvec, decltype(xx.row(i))>(xx.row(i), (size_t)ss[i], replace[i], parallel);
            }
        }
        return res;
    }

    template <class RngFunction, class... Args>
    Rcpp::NumericMatrix colRng(size_t nrow, size_t ncol, RngFunction rngF, Args &&...args)
    {
        Rcpp::NumericMatrix res(nrow, ncol);
        arma::mat Res = rangen_internal::getArmaFrom(res, false);

        if constexpr (sizeof...(Args) == 0)
        {
            for (size_t i = 0; i < ncol; ++i)
            {
                Res.col(i) = rngF(nrow);
            }
        }
        else
        {
            for (size_t i = 0; i < ncol; ++i)
            {
                Res.col(i) = rngF(nrow, (args[i])...); // Unpack i-th element of each arg
            }
        }
        return res;
    }

}