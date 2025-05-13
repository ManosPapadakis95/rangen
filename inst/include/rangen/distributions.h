#pragma once

#include <RcppArmadillo.h>
#include "Random.h"
#include "assertions.hpp"

namespace rangen
{

	namespace rangen_internal
	{

		template <typename T>
		inline constexpr bool is_arma = std::is_base_of<arma::Mat<typename T::elem_type>, T>::value ||
										std::is_base_of<arma::Col<typename T::elem_type>, T>::value ||
										std::is_base_of<arma::Row<typename T::elem_type>, T>::value;

		template <typename T>
		inline constexpr bool is_rcpp = std::is_same<Rcpp::NumericMatrix, T>::value ||
										std::is_same<Rcpp::IntegerMatrix, T>::value ||
										std::is_same<Rcpp::CharacterMatrix, T>::value ||
										std::is_same<Rcpp::StringMatrix, T>::value ||
										std::is_same<Rcpp::NumericVector, T>::value ||
										std::is_same<Rcpp::IntegerVector, T>::value ||
										std::is_same<Rcpp::CharacterVector, T>::value ||
										std::is_same<Rcpp::StringVector, T>::value;

		template <class T>
		size_t nrow(T x)
		{
			Assertion::has_value_type<T>::check_concept();

			size_t res;

			if constexpr (is_arma<T>)
			{
				res = x.n_rows;
			}
			else if constexpr (is_rcpp<T>)
			{
				res = x.nrow();
			}
			return res;
		}

		template <class T>
		size_t ncol(T x)
		{
			Assertion::has_value_type<T>::check_concept();

			size_t res;

			if constexpr (std::is_base_of<arma::Mat<typename T::value_type>, T>::value)
			{
				res = x.n_cols;
			}
			else if constexpr (is_rcpp<T>)
			{
				res = x.ncol();
			}
			return res;
		}

		template <class T>
		T getVector(size_t size)
		{

			Assertion::has_subscript_operator<T>::check_concept();
			Assertion::has_value_type<T>::check_concept();
			Assertion::has_size<T>::check_concept();

			T res;

			if constexpr (is_arma<T>)
			{
				res = T(size, arma::fill::none);
			}
			else if constexpr (is_rcpp<T>)
			{
				res = T(size);
			}
			return res;
		}

		template <class T>
		T getMatrix(size_t nrow, size_t ncol)
		{

			Assertion::has_subscript_operator<T>::check_concept();
			Assertion::has_value_type<T>::check_concept();
			Assertion::has_size<T>::check_concept();

			T res;

			if constexpr (is_arma<T>)
			{
				res = T(nrow, ncol, arma::fill::none);
			}
			else if constexpr (is_rcpp<T>)
			{
				res = T(nrow, ncol);
			}
			return res;
		}

		template <class T>
		arma::Mat<typename std::remove_reference<typename T::value_type>::type> getArmaFrom(T x, const bool copy = true)
		{

			Assertion::has_subscript_operator<T>::check_concept();
			Assertion::has_size<T>::check_concept();

			arma::Mat<typename std::remove_reference<typename T::value_type>::type> res(x.begin(), nrow(x), ncol(x), copy);
			return res;
		}

		template <class T, class Generator, class... Args>
		T generic(size_t size, Args... args)
		{
			T res = rangen_internal::getVector<T>(size);

			Generator rng(args...);
			for (size_t i = 0; i < size; ++i)
			{
				res[i] = rng();
			}
			return res;
		}
	}

	template <class T>
	T sample(T x, size_t size, const bool replace = false)
	{
		T res = rangen_internal::getVector<T>(size);

		if constexpr(std::is_integral_v<typename T::value_type>){
			if(x.size() == 1){
				const size_t n = x[0];
				for (unsigned int i = 0; i < n; ++i){
					res[i] = i + 1;
				}
			}
		}

		if (replace)
		{
			uniform<integer, true> rng(0, size - 1);
			for (unsigned int i = 0; i < size; ++i)
			{
				res[i] = x[rng()];
			}
		}
		else
		{
			for (unsigned int i = 0; i < size; ++i)
			{
				res[i] = x[rng()];
			}
		}
		return res;
	}
	
	// template<class T, class L>
	// T colSample(T x, size_t size, L replace){
	// 	Assertion::has_value_type<T>::check_concept();
	// 	Assertion::has_value_type<T>::check_concept();
	// 	Assertion::has_col<T>::check_concept();
	// 	Assertion::has_col<L>::check_concept();

	// 	size_t n = ncol(x);
	// 	T res = rangen_internal::getVector<T>(n);

	// 	for(size_t i=0; i<n; ++i){
	// 		if constexpr (is_arma<T>) {
	// 			res.col(i) = sample<arma::Col<typename T::vale_type>>(x.col(i), size, replace[i]);
	// 		} else if constexpr (is_rcpp<T> || is_custom<T>) {
	// 			res.column(i) = sample>(x.column(i), size, replace[i]);
	// 		}
	// 	}
	// }

	template <class L>
	Rcpp::NumericMatrix colSample(Rcpp::NumericMatrix x, size_t size, L replace)
	{
		size_t n = rangen_internal::ncol(x);
		arma::mat xx = rangen_internal::getArmaFrom(x, false);
		Rcpp::NumericMatrix res(xx.n_rows, xx.n_cols);
		arma::mat Res = rangen_internal::getArmaFrom(res, false);

		for (size_t i = 0; i < n; ++i)
		{
			Res.col(i) = sample<arma::colvec>(xx.col(i), size, replace[i]);
		}
		return res;
	}

	template <class L>
	Rcpp::NumericMatrix rowSample(Rcpp::NumericMatrix x, size_t size, L replace)
	{
		size_t n = rangen_internal::nrow(x);
		arma::mat xx = rangen_internal::getArmaFrom(x, false);
		Rcpp::NumericMatrix res(xx.n_rows, xx.n_cols);
		arma::mat Res = rangen_internal::getArmaFrom(res, false);

		for (size_t i = 0; i < n; ++i)
		{
			Res.row(i) = sample<arma::rowvec>(xx.row(i), size, replace[i]);
		}
		return res;
	}

	template <class T>
	T runif(size_t size, double min = 0.0, double max = 1.0)
	{
		return rangen_internal::generic<T, uniform<real>>(size, min, max);
	}

	template <class T>
	T rbeta(size_t size, double alpha, double beta)
	{
		return rangen_internal::generic<T, Beta>(size, alpha, beta);
	}

	template <class T>
	T rexp(size_t size, double rate = 1.0)
	{
		return rangen_internal::generic<T, Exp>(size, rate);
	}

	template <class T>
	T rchisq(size_t size, double df)
	{
		return rangen_internal::generic<T, Chisq>(size, df);
	}

	template <class T>
	T rgamma(size_t size, double shape, double rate = 1.0)
	{
		return rangen_internal::generic<T, Gamma>(size, shape, rate);
	}

	template <class T>
	T rgeom(size_t size, double prob)
	{
		return rangen_internal::generic<T, Geom>(size, prob);
	}

	template <class T>
	T rcauchy(size_t n, double location = 0.0, double scale = 1.0)
	{
		return rangen_internal::generic<T, Cauchy>(n, location, scale);
	}

	template <class T>
	T rt(size_t n, double df, double ncp)
	{
		return rangen_internal::generic<T, StudentT>(n, df, ncp);
	}

	template <class T>
	T rpareto(size_t n, double shape = 1.0, double scale = 1.0)
	{
		return rangen_internal::generic<T, Pareto>(n, shape, scale);
	}

	template <class T>
	T rfrechet(size_t n, double shape = 1.0, double mean = 0.0, double scale = 1.0)
	{
		return rangen_internal::generic<T, Frechet>(n, shape, mean, scale);
	}

	template <class T>
	T rlaplace(size_t n, double mean = 0.0, double sigma = 1.0)
	{
		return rangen_internal::generic<T, Laplace>(n, mean, sigma);
	}

	template <class T>
	T rgumble(size_t n, double mean = 0.0, double sigma = 1.0)
	{
		return rangen_internal::generic<T, Gumbel>(n, mean, sigma);
	}

	template <class T>
	T rarcsine(size_t n, double min = 0.0, double max = 1.0)
	{
		return rangen_internal::generic<T, Arcsine>(n, min, max);
	}

}