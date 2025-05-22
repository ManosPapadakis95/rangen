#pragma once

#include "Random.h"
#include "assertions.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rangen
{

	namespace rangen_internal
	{
		inline unsigned int get_num_of_threads()
		{
#ifdef _OPENMP
			return omp_get_max_threads();
#else
			return 0;
#endif
		}

		template <typename T>
		inline constexpr bool is_arma_mat =
#ifdef ARMA_VERSION_MAJOR
			std::is_base_of<arma::Mat<typename T::elem_type>, T>::value ||
			std::is_base_of<arma::subview<typename T::elem_type>, T>::value
#else
			false
#endif
			;

		template <typename T>
		inline constexpr bool is_arma_vec =
#ifdef ARMA_VERSION_MAJOR
			std::is_base_of<arma::Col<typename T::elem_type>, T>::value ||
			std::is_base_of<arma::Row<typename T::elem_type>, T>::value ||
			std::is_base_of<arma::subview_col<typename T::elem_type>, T>::value ||
			std::is_base_of<arma::subview_row<typename T::elem_type>, T>::value
#else
			false
#endif
			;

		template <typename T>
		inline constexpr bool is_arma = is_arma_mat<T> || is_arma_vec<T>;

		template <typename T>
		inline constexpr bool is_rcpp_mat =
#ifdef ARMA_VERSION_MAJOR
			std::is_same<Rcpp::NumericMatrix, T>::value ||
			std::is_same<Rcpp::IntegerMatrix, T>::value ||
			std::is_same<Rcpp::CharacterMatrix, T>::value ||
			std::is_same<Rcpp::StringMatrix, T>::value
#else
			false
#endif
			;

		template <typename T>
		inline constexpr bool is_rcpp_vec =
#ifdef ARMA_VERSION_MAJOR
			std::is_same<Rcpp::NumericVector, T>::value ||
			std::is_same<Rcpp::IntegerVector, T>::value ||
			std::is_same<Rcpp::CharacterVector, T>::value ||
			std::is_same<Rcpp::StringVector, T>::value
#else
			false
#endif
			;

		template <typename T>
		inline constexpr bool is_rcpp = is_rcpp_mat<T> || is_rcpp_vec<T>;

		template <class T>
		size_t nrow(T x)
		{
			size_t res;

			if constexpr (is_arma_mat<T>)
			{
				res = x.n_rows;
			}
			else if constexpr (is_rcpp_mat<T>)
			{
				res = x.nrow();
			}
			else if constexpr (!is_arma_mat<T> && !is_rcpp_mat<T>)
			{
				res = x.nrow;
			}
			return res;
		}

		template <class T>
		size_t ncol(T x)
		{
			size_t res;

			if constexpr (is_arma_mat<T>)
			{
				res = x.n_cols;
			}
			else if constexpr (is_rcpp_mat<T>)
			{
				res = x.ncol();
			}
			else if constexpr (!is_arma_mat<T> && !is_rcpp_mat<T>)
			{
				res = x.ncol;
			}
			return res;
		}

		template <class T>
		size_t size(T x)
		{
			size_t res;

			if constexpr (is_arma_vec<T>)
			{
				res = x.n_elem;
			}
			else if constexpr (is_rcpp_vec<T>)
			{
				res = x.size();
			}
			else if constexpr (!is_arma_vec<T> && !is_rcpp_vec<T>)
			{
				res = x.size;
			}
			return res;
		}

		template <class T>
		T getVector(size_t size)
		{

			Assertion::has_subscript_operator<T>::check_concept();
			Assertion::has_size<T>::check_concept();

			T res;

			if constexpr (is_arma_vec<T>)
			{
				res = T(size, arma::fill::none);
			}
			else if constexpr (is_rcpp_vec<T>)
			{
				res = T(size);
			}
			else if constexpr (!is_arma_vec<T> && !is_rcpp_vec<T>)
			{
				res = T(size);
			}
			return res;
		}

		template <class T>
		T getMatrix(size_t nrow, size_t ncol)
		{

			Assertion::has_subscript_operator<T>::check_concept();
			Assertion::has_size<T>::check_concept();

			T res;

			if constexpr (is_arma_mat<T>)
			{
				res = T(nrow, ncol, arma::fill::none);
			}
			else if constexpr (is_rcpp_mat<T>)
			{
				res = T(nrow, ncol);
			}
			else if constexpr (!is_arma_mat<T> && !is_rcpp_mat<T>)
			{
				res = T(nrow, ncol);
			}
			return res;
		}

		template <class T, class Generator>
		T generic(size_t size, Generator rng)
		{
			T res = rangen_internal::getVector<T>(size);

			for (size_t i = 0; i < size; ++i)
			{
				res[i] = rng();
			}
			return res;
		}

		template <class T, class Generator, class... Args>
		T generic(size_t size, Args... args)
		{
			Generator rng(args...);
			return generic<T>(size, rng);
		}
	}

	template <class Ret, class T = Ret>
	Ret sample(T x, size_t size, const bool replace = false, const bool thread_safe = false)
	{
		Ret res = rangen_internal::getVector<Ret>(size);
		size_t n = rangen_internal::size(x);

		if (replace)
		{
			if(thread_safe){
				uniform<integer, true> rng(0, n - 1);
				for (size_t i = 0; i < size; ++i)
				{
					res[i] = x[rng()];
				}
			}else{
				irng_rep.set_bounds(0,n-1);
				for (size_t i = 0; i < size; ++i)
				{
					res[i] = x[irng_rep()];
				}
			}
		}
		else
		{
			if(thread_safe){
				uniform<integer> rng(0, n - 1);
				for (size_t i = 0; i < size; ++i)
				{
					res[i] = x[rng()];
				}
			}else{
				irng.set_bounds(0,n-1);
				for (size_t i = 0; i < size; ++i)
				{
					res[i] = x[irng()];
				}
			}
		}
		return res;
	}

	template <class Ret>
	Ret sample(size_t x, size_t size, const bool replace = false, const bool thread_safe = false)
	{
		Ret res = rangen_internal::getVector<Ret>(size);
		for (size_t i = 0; i < x; ++i)
		{
			res[i] = i + 1;
		}
		return sample<Ret>(res, size, replace, thread_safe);
	}

	template <class Ret, class T, class S, class L>
	Ret colSample(T x, S size, L replace, const bool parallel = false, const size_t cores = rangen_internal::get_num_of_threads())
	{
		const size_t n = rangen_internal::ncol(x);
		const size_t m = *std::max_element(size.begin(), size.end());
		Ret res = rangen_internal::getMatrix<Ret>(m, n);

		if (parallel)
		{

			if constexpr (rangen_internal::is_rcpp_mat<T>)
			{
				Rcpp::stop("Rcpp matrices are not thread safe. Please use `parallel = false`.");
			}

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (size_t i = 0; i < n; ++i)
			{
				if constexpr (rangen_internal::is_arma_mat<T>)
				{
					res.col(i) = sample<arma::Col<typename T::value_type>>(x.col(i), size[i], replace[i], parallel);
				}
				else if constexpr (rangen_internal::is_rcpp_mat<T>)
				{
					res.column(i) = sample<decltype(x.column(i))>(x.column(i), size[i], replace[i], parallel);
				}
				else if constexpr (!rangen_internal::is_arma_mat<T> && !rangen_internal::is_rcpp_mat<T>)
				{
					res.column(i) = sample<typename T::col_type>(x.column(i), size[i], replace[i], parallel);
				}
			}
		}
		else
		{
			for (size_t i = 0; i < n; ++i)
			{
				if constexpr (rangen_internal::is_arma_mat<T>)
				{
					res.col(i) = sample<arma::Col<typename T::value_type>>(x.col(i), size[i], replace[i], parallel);
				}
				else if constexpr (rangen_internal::is_rcpp_mat<T>)
				{
					res.column(i) = sample<decltype(x.column(i))>(x.column(i), size[i], replace[i], parallel);
				}
				else if constexpr (!rangen_internal::is_arma_mat<T> && !rangen_internal::is_rcpp_mat<T>)
				{
					res.column(i) = sample<typename T::col_type>(x.column(i), size[i], replace[i], parallel);
				}
			}
		}
		return res;
	}

	template <class Ret, class T, class S, class L>
	Ret rowSample(T x, S size, L replace, const bool parallel = false, const size_t cores = rangen_internal::get_num_of_threads())
	{
		const size_t m = rangen_internal::nrow(x);
		const size_t n = *std::max_element(size.begin(), size.end());
		Ret res = rangen_internal::getMatrix<Ret>(m, n);

		if (parallel)
		{

			if (rangen_internal::is_rcpp_mat<T>)
			{
				Rcpp::stop("Rcpp matrices are not thread safe. Please use `parallel = false`.");
			}

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (size_t i = 0; i < m; ++i)
			{
				if constexpr (rangen_internal::is_arma_mat<T>)
				{
					res.row(i) = sample<arma::Row<typename T::value_type>>(x.row(i), size[i], replace[i], parallel);
				}
				else if constexpr (rangen_internal::is_rcpp_mat<T>)
				{
					res.row(i) = sample<decltype(x.row(i))>(x.row(i), size[i], replace[i], parallel);
				}
				else if constexpr (!rangen_internal::is_arma_mat<T> && !rangen_internal::is_rcpp_mat<T>)
				{
					res.row(i) = sample<typename T::row_type>(x.row(i), size[i], replace[i], parallel);
				}
			}
		}
		else
		{
			for (size_t i = 0; i < m; ++i)
			{
				if constexpr (rangen_internal::is_arma_mat<T>)
				{
					res.row(i) = sample<arma::Row<typename T::value_type>>(x.row(i), size[i], replace[i], parallel);
				}
				else if constexpr (rangen_internal::is_rcpp_mat<T>)
				{
					res.row(i) = sample<decltype(x.row(i))>(x.row(i), size[i], replace[i], parallel);
				}
				else if constexpr (!rangen_internal::is_arma_mat<T> && !rangen_internal::is_rcpp_mat<T>)
				{
					res.row(i) = sample<typename T::row_type>(x.row(i), size[i], replace[i], parallel);
				}
			}
		}
		return res;
	}

	template <class T>
	T runif(size_t size, double min = 0.0, double max = 1.0)
	{
		rng2.set_bounds(min, max);
		return rangen_internal::generic<T>(size, rng2);
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