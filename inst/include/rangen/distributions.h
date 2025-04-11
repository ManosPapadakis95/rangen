#pragma once


#include <RcppArmadillo.h>
#include "Random.h"
#include "assertions.hpp"

namespace rangen
{

	namespace rangen_internal
	{
		template<class T>
		T getOptimType(size_t size){

			Assertion::has_subscript_operator<T>::check_concept();
			Assertion::has_value_type<T>::check_concept();
			Assertion::has_size<T>::check_concept();
			
			T res;
			
			if constexpr(std::is_base_of<arma::Mat<typename T::value_type>, T>::value){
				res = T(size, arma::fill::none);
			}else{
				res = T(size);
			}
			return res;
		}

		template<typename T>
		inline constexpr bool is_arma = std::is_base_of<arma::Mat<typename T::value_type>, T>::value;

		template<typename T>
		inline constexpr bool is_rcpp = std::is_same<Rcpp::NumericMatrix, T>::value || 
		std::is_same<Rcpp::IntegerMatrix, T>::value || 
		std::is_same<Rcpp::CharacterMatrix, T>::value || 
		std::is_same<Rcpp::StringMatrix, T>::value;

		template<typename T>
		inline constexpr bool is_custom = Assertion::has_nrow<T>::check_concept();

		template<class T>
		size_t nrow(T x) {
			Assertion::has_value_type<T>::check_concept();

			if constexpr (is_arma<T>) {
				return x.n_rows;
			} else if constexpr (is_rcpp<T> || is_custom<T>) {
				return x.nrow();
			}
		}

		
		template<class T>
		size_t ncol(T x) {
			Assertion::has_value_type<T>::check_concept();

			if constexpr (std::is_base_of<arma::Mat<typename T::value_type>, T>::value) {
				return x.n_cols;
			} else if constexpr (
				std::is_same<Rcpp::NumericMatrix, T>::value || 
				std::is_same<Rcpp::IntegerMatrix, T>::value || 
				std::is_same<Rcpp::CharacterMatrix, T>::value || 
				std::is_same<Rcpp::StringMatrix, T>::value ||
				Assertion::has_ncol<T>::check_concept()) {
				return x.ncol();
			}
		}


		template<class T, class Generator, class ...Args>
		T generic(size_t size, Args... args)
		{
			T res = rangen_internal::getOptimType<T>(size);

			Generator rng(args...);
			for (size_t i = 0; i < size; ++i)
			{
				res[i] = rng();
			}
			return res;
		}
	}

	template<class T>
	T sample(size_t n, size_t size, const bool replace = false)
	{
        Assertion::has_value_int<T>();
     
        T res = rangen_internal::getOptimType<T>(size);

		if (replace)
		{
			uniform<integer, true> rng(1, n);
			for (unsigned int i = 0; i < size; ++i)
			{
				res[i] = rng();
			}
		}
		else
		{
			uniform<integer> rng(1, n);
			for (unsigned int i = 0; i < size; ++i)
			{
				res[i] = rng();
			}
		}
		return res;
	}

	template<class T>
	T sample(T x, size_t size, const bool replace = false)
	{
		T res = rangen_internal::getOptimType<T>(size);
		
		if (replace)
		{
			uniform<integer, true> rng(0, x.size() - 1);
			for (unsigned int i = 0; i < size; ++i)
			{
				res[i] = x[rng()];
			}
		}
		else
		{
			uniform<integer> rng(0, x.size() - 1);
			for (unsigned int i = 0; i < size; ++i)
			{
				res[i] = x[rng()];
			}
		}
		return res;
	}

    template<class T, class L>
	colSample(T x, size_t size, L replace){
		Assertion::has_value_type<T>::check_concept();
		Assertion::has_value_type<T>::check_concept();
		Assertion::has_col<T>::check_concept();
		Assertion::has_col<L>::check_concept();

		size_t n = ncol(x);
		T res = rangen_internal::getOptimType<T>(n);

		for(size_t i=0; i<n; ++i){
			if constexpr (is_arma<T>) {
				x.col(i) = sample<arma::Col<typename T::vale_type>>(x.col(i), size, replace[i]);
			} else if constexpr (is_rcpp<T> || is_custom<T>) {
				x.column(i) = sample>(x.column(i), size, replace[i]);
			}
		}
		
		return res;
	}

    template<class T>
	T runif(size_t size, const double min, const double max)
	{
		return rangen_internal::<T,uniform<real>>(size,min,max);
	}

	template<class T>
	T rbeta(size_t size, double alpha, double beta)
	{
		return rangen_internal::<T,Beta>(size,alpha,beta);
	}
	
	template<class T>
	T rexp(size_t size, double rate)
	{
		return rangen_internal::<T,Exp>(size,rate);
	}
	
	template<class T>
	T rchisq(size_t size, double df)
	{
		return rangen_internal::<T,Chisq>(size,df);
	}
	
	template<class T>
	T rgamma(size_t size, double shape, double rate)
	{
		return rangen_internal::<T,Gamma>(size,shape, rate);
	}
	
	template<class T>
	T rgeom(size_t size, double prob)
	{
		return rangen_internal::<T,Geom>(size,prob);
	}
	
	template<class T>
	T rcauchy(size_t size, double location, double scale)
	{
		return rangen_internal::<T,Cauchy>(size,location,scale);
	}
	
	template<class T>
	T rt(size_t size, double df, double ncp)
	{
		return rangen_internal::<T,StudentT>(size,df, ncp);
	}

	template <class T>
	T rfrechet(size_t size, double shape, double mean, double scale) {
		return rangen_internal::<T, Frechet>(size, shape, mean, scale);
	}

	template <class T>
	T rlaplace(size_t size, double mean, double sigma) {
		return rangen_internal::<T, Laplace>(size, mean, sigma);
	}

	template <class T>
	T rgumble(size_t size, double mean, double sigma) {
		return rangen_internal::<T, Gumbel>(size, mean, sigma);
	}

	template <class T>
	T rarcsine(size_t size, double min, double max) {
		return rangen_internal::<T, Arcsine>(size, min, max);
	}





}