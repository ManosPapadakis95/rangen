#pragma once


#include <RcppArmadillo.h>
#include "Random.h"
#include "assertions.hpp"

namespace ragen
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

	template<class T>
	T sample(size_t n, size_t size, const bool replace)
	{
        Assertion::has_value_int<T>();
     
        T res = getOptimType<T>(size);

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
	T sample(T x, size_t size, const bool replace)
	{
		T res = getOptimType<T>(size);
		
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
    
    template<class T, class Generator, class ...Args>
	T generic(size_t size, Args... args)
	{
		T res = getOptimType<T>(size);

		Generator rng(args...);
		for (size_t i = 0; i < size; ++i)
		{
			res[i] = rng();
		}
		return res;
	}

    template<class T>
	T runif(size_t size, const double min, const double max)
	{
		return generic<T,uniform<real>>(size,min,max);
	}

	template<class T>
	T rbeta(size_t size, double alpha, double beta)
	{
		return generic<T,Beta>(size,alpha,beta);
	}
	
	template<class T>
	T rexp(size_t size, double rate)
	{
		return generic<T,Exp>(size,rate);
	}
	
	template<class T>
	T rchisq(size_t size, double df)
	{
		return generic<T,Chisq>(size,df);
	}
	
	template<class T>
	T rgamma(size_t size, double shape, double rate)
	{
		return generic<T,Gamma>(size,shape, rate);
	}
	
	template<class T>
	T rgeom(size_t size, double prob)
	{
		return generic<T,Geom>(size,prob);
	}
	
	template<class T>
	T rcauchy(size_t size, double location, double scale)
	{
		return generic<T,Cauchy>(size,location,scale);
	}
	
	template<class T>
	T rt(size_t size, double df, double ncp)
	{
		return generic<T,StudentT>(size,df, ncp);
	}

}