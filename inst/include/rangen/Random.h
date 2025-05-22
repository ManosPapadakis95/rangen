#ifndef RANDOM_GENERATORS_H

#define RANDOM_GENERATORS_H
#include <cmath>
#include <limits>
#include <type_traits>
#include <chrono>
#include <vector>
#include <algorithm>

//[[Rcpp::depends(zigg, RcppArmadillo)]]

#include <zigg/header>

namespace rangen
{

	using std::numeric_limits;
	using real = std::false_type;
	using integer = std::true_type;
	using std::iota;
	using std::vector;
	static zigg::Ziggurat ziggurat;

	namespace internal
	{

		inline long long int get_cur_nano()
		{
			return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		}

		class Integer_Core
		{
		public:
			// The next 3 lines are the requirements for UniformRandomBitGenerator.
			using result_type = uint32_t;
			static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }
			static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

			using type = uint64_t;

			inline void setSeed(type state){
				rng.state = state;
			}

		protected:
			struct pcg32_random_t
			{
				type state;
				type inc;
				pcg32_random_t(type init = get_cur_nano()) : state(init), inc(init) {}
			};

			pcg32_random_t rng;

			Integer_Core(type init = get_cur_nano()) : rng(pcg32_random_t(init)) {}

			result_type pcg32_random_r()
			{
				type oldstate = rng.state;
				// Advance internal state
				rng.state = oldstate * 6364136223846793005ULL + (rng.inc | 1);
				// Calculate output function (XSH RR), uses old state for max ILP
				result_type xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
				result_type rot = oldstate >> 59u;
				return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
			}
		};
	}

	template <typename T, bool replace = false>
	class uniform : public internal::Integer_Core
	{

		vector<size_t> indices;

		void remove_index(result_type i)
		{
			this->indices[(size_t)i] = this->indices.back();
			this->indices.pop_back();
		}

	public:

		uniform(int32_t min_bound = 0, int32_t max_bound = 0)
		{
			this->set_bounds(min_bound, max_bound);
			this->indices.resize(std::abs(max_bound - min_bound + 1));
			iota(this->indices.begin(), this->indices.end(), min_bound);
		}

		void set_bounds(int32_t min_bound, int32_t max_bound){
			this->indices.resize(std::abs(max_bound - min_bound + 1));
			iota(this->indices.begin(), this->indices.end(), min_bound);
		}

		result_type operator()()
		{
			auto index = (size_t)this->pcg32_random_r() % this->indices.size();
			auto res = this->indices[index];
			this->remove_index(index);
			return res;
		}
	};

	template <typename T>
	class uniform<T, true> : public internal::Integer_Core
	{
		// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
		// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

		int32_t min_bound, max_bound;

	public:
		uniform(int32_t min_bound = 0, int32_t max_bound = 0) : min_bound(min_bound), max_bound(max_bound) {}

		void set_bounds(int32_t min_bound, int32_t max_bound){
			this->min_bound = min_bound; 
			this->max_bound = max_bound;
		}

		int32_t operator()()
		{
			return this->pcg32_random_r() % this->max_bound + this->min_bound;
		}
	};

	template <>
	class uniform<real, false> : public internal::Integer_Core
	{

		double min, max;

	public:
		uniform(double min = 0.0, double max = 1.0) : min(min), max(max) {}

		void set_bounds(double min, double max){
			this->min = min; 
			this->max = max;
		}

		inline double operator()()
		{
			return min + (this->pcg32_random_r() * (max - min) / internal::Integer_Core::max());
		}
	};

	inline uniform<integer> irng;
	inline uniform<integer, true> irng_rep;
	inline uniform<real> rng(0, 1);
	inline uniform<real> rng2;

	class Gamma
	{
		double init_shape, shape, rate, d, c, boosted_shape, inv_init_shape;
		bool boosted;

	public:
		Gamma(double shape, double rate = 1.0) : shape(shape), rate(1.0 / rate)
		{
			if (shape < 1.0)
			{
				boosted_shape = shape + 1.0;
				boosted = true;
				inv_init_shape = 1.0 / shape; // Precompute to avoid division
			}
			else
			{
				boosted_shape = shape;
				boosted = false;
			}
			d = boosted_shape - 1.0 / 3.0;
			c = 1.0 / std::sqrt(9.0 * d);
		}

		double operator()()
		{
			while (true)
			{
				double x = ziggurat.rnorm();
				double v = 1.0 + c * x;
				v = v * v * v;
				if (v > 0)
				{
					double u = rng();
					double x2 = x * x;
					if (u < 1.0 - 0.0331 * x2 * x2 || std::log(u) < 0.5 * x2 + d * (1.0 - v + std::log(v)))
					{
						double res = d * v * rate;
						if (boosted)
						{
							return res * std::exp(std::log(u) * inv_init_shape); // Faster than pow()
						}
						else
						{
							return res;
						}
					}
				}
			}
		}
	};

	class BetaOne : public Gamma
	{

	public:
		BetaOne(double alpha) : Gamma(alpha, 1.0) {}

		inline double operator()()
		{
			double x = Gamma::operator()();
			return x / (x + Gamma::operator()());
		}
	};

	class Beta : public BetaOne
	{
		double beta;
		Gamma beta_d;

	public:
		Beta(double alpha, double beta) : BetaOne(alpha), beta(beta), beta_d(Gamma(beta, 1.0)) {}

		inline double operator()()
		{
			double x = Gamma::operator()();
			return x / (x + beta_d());
		}
	};

	class Exp : public Gamma
	{
	public:
		Exp(double rate = 1.0) : Gamma(1.0, rate) {}
	};

	class Chisq : public Gamma
	{
	public:
		// rate must be 2 but Gamma divides with 1. So to undo it we need to divided first with 1 and pass it.
		Chisq(double df) : Gamma(df * 0.5, 0.5) {}
	};

	class Geom
	{
		double lambda;

	public:
		Geom(double prob) : lambda(std::log(1 - prob)) {}

		inline double operator()()
		{
			return std::floor(std::log(rng()) / lambda);
		}
	};

	class Cauchy
	{
		double location, scale;

	public:
		Cauchy(double location = 0.0, double scale = 1.0) : location(location), scale(scale) {}

		double operator()()
		{
			return location + scale * std::tan(M_PI * ziggurat.rnorm());
		}
	};

	class StudentT : public Chisq
	{
		double sqrt_df, ncp_sqrt_df;

	public:
		StudentT(double df, double ncp) : Chisq(df), sqrt_df(std::sqrt(df)), ncp_sqrt_df(ncp * sqrt_df) {}

		double operator()()
		{
			// return (ziggurat.rnorm() + ncp) / std::sqrt(Chisq::operator()() / df);
			return (ziggurat.rnorm() * sqrt_df + ncp_sqrt_df) / std::sqrt(Chisq::operator()());
		}
	};
	
	class Pareto
	{
		double shape, scale;

	public:
		Pareto(double shape = 1.0, double scale = 1.0) : shape(shape), scale(scale) {}

		double operator()()
		{
			return scale*std::pow(1-rng(), -1/shape);
		}
	};

	class Frechet {
		double shape, mean, scale;
	
	public:
		Frechet(double shape = 1.0, double mean = 0.0, double scale = 1.0) : shape(shape), mean(mean), scale(scale) {}
	
		double operator()() {
			return mean + scale * std::pow(-std::log(rng()), 1.0 / shape);
		}
	};
	
	class Laplace {
		double mean, sigma;
	
	public:
		Laplace(double mean = 0.0, double sigma = 1.0) : mean(mean), sigma(sigma) {}
	
		double operator()() {
			double u = rng() - 0.5;
			return mean - sigma * ((u > 0) - (u < 0)) * std::log(1 - 2 * std::abs(u));
		}
	};
	
	class Gumbel {
		double mean, sigma;
	
	public:
		Gumbel(double mean = 0.0, double sigma =  1.0) : mean(mean), sigma(sigma) {}
	
		double operator()() {
			return mean - sigma * std::log(-std::log(rng()));
		}
	};
	
	class Arcsine {
		double min, max;
	
	public:
		Arcsine(double min = 0.0, double max = 1.0) : min(min), max(max) {}
	
		double operator()() {
			constexpr double P = M_PI / 2.0;
			double k = std::sin(rng() * P);
			return min + (max - min) * k * k;
		}
	};

	inline void setSeed(const size_t s){
		ziggurat.setSeed(static_cast<uint32_t>(s));
		rng.setSeed(static_cast<uint64_t>(s));
		rng2.setSeed(static_cast<uint64_t>(s));
		irng.setSeed(static_cast<uint64_t>(s));
		irng_rep.setSeed(static_cast<uint64_t>(s));
	}
	
}

#endif