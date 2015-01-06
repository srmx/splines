#include <stdexcept>
#include <assert.h>
#include <iostream>
#include <algorithm>

#include "splines.h"


/*
 * k!
 */
constexpr int fact(int k)
{
	return k > 0 ? k * fact(k-1) : 1;
}

/*
 * n! / k! (n-k)! for k=0..n
 */
constexpr int bin_coef(int n, int k)
{
	return k >= 0 && k <= n ? fact(n) / (fact(k) * fact(n-k)) : 0;
}

/*
 * non-uniform knots spline
 */
template <int n /* degree of the spline */>
class non_uniform_spline : 
	public basic_spline
{
private:
	std::vector<double> knots;
	std::vector<double> coefs;

private:
	/*
	 * $ \sum N_{i}^{n}C_{i}=A_{1,0}\cdot\left\{ A_{2,0}\cdot\left\{ A_{3,0}\cdot\left\{ ...\right\} +B_{3,0}\cdot\left\{ ...\right\} \right\} +B_{2,0}\cdot\left\{ A_{3,1}\cdot\left\{ ...\right\} +B_{3,1}\cdot\left\{ ...\right\} \right\} \right\} +  B_{1,0}\cdot\left\{ A_{2,1}\cdot\left\{ A_{3,1}\cdot\left\{ ...\right\} +B_{3,1}\cdot\left\{ ...\right\} \right\} +B_{2,1}\cdot\left\{ A_{3,2}\cdot\left\{ ...\right\} +B_{3,2}\cdot\left\{ ...\right\} \right\} \right\} $
	 * where
	 * $ A_{k,j}=\frac{x-u_{l-j}}{u_{l+k-j}-u_{l-j}} $
	 * $ B_{k,j}=\frac{u_{l-j+k}-x}{u_{l-j+k}-u_{l-j}} $, $k = 1..n, j = 0..n-1$
	 */
	template <int k, int j>
	inline double sum(int l, double x) const
	{
		if (k == n + 1)
			return coefs[l - j];

		constexpr int nk = k > n ? k : k + 1;
		constexpr int nj = k > n ? j : j + 1;

		double const Akj = (x - knots[l - j]) / (knots[l - j + k] - knots[l - j]);
		double const Bkj = 1.0 - Akj;

		return Akj * sum<nk,j>(l, x) + Bkj * sum<nk,nj>(l, x);	
	}

	inline double sum(int l, double x) const
	{
		return sum<1,0>(l, x);
	}

	/*
	 * $\frac{d}{dx}N_{i}^{n}=Q_{0,i}N_{i}^{n-1}+Q_{1,i}N_{i+1}^{n-1}$
	 * $Q_{0,i}=\frac{n}{u_{i+n}-u_{i}}$
	 * $Q_{1,i}=-\frac{n}{u_{i+n+1}-u_{i+1}}$
	 * for $k=n$ we return 
	 * $Q_{0,i}C_{i}+Q_{1,i-1}C_{i-1}$
 	 */
	template <int k, int j>
	inline double der_sum(int l, double x) const
	{
		if (k == n)
			return (coefs[l - j] - coefs[l - j - 1]) / (knots[l - j + k] - knots[l - j]);

		constexpr int nk = k > n ? k : k + 1;
		constexpr int nj = k > n ? j : j + 1;

		double const Akj = (x - knots[l - j]) / (knots[l - j + k] - knots[l - j]);
		double const Bkj = 1.0 - Akj;

		return Akj * der_sum<nk,j>(l, x) + Bkj * der_sum<nk,nj>(l, x);
	}

	inline double der_sum(int l, double x) const
	{
		assert(n >= 1);
		return der_sum<1,0>(l, x) * n;
	}

	/*
	 * for k = n-1:
	 * $\frac{d^{2}}{dx^{2}}N_{i}^{n}=Q_{0,i}N_{i}^{n-2}+Q_{1,i}N_{i+1}^{n-2}+Q_{2,i}N_{i+2}^{n-2}$
	 * $Q_{0,i}=\frac{n}{u_{i+n}-u_{i}}\cdot\frac{n-1}{u_{i+n-1}-u_{i}}$
	 * $Q_{1,i}=-\frac{n}{u_{i+n}-u_{i+1}}\cdot\frac{n-1}{u_{i+n}-u_{i}}-\frac{n}{u_{i+n}-u_{i+1}}\cdot\frac{n-1}{u_{i+n+1}-u_{i+1}}$
	 * $Q_{2,i}=\frac{n}{u_{i+n+1}-u_{i+1}}\cdot\frac{n-1}{u_{i+n+1}-u_{i+2}}$
	 */
	template <int k, int j>
	inline double der2_sum(int l, double x) const
	{
		if (k == n - 1)
		{
			int const i = l - j;

			double const Q0 = 1.0 / ((knots[i+n] - knots[i]) * (knots[i+n-1] - knots[i]));
			double const Q2 = 1.0 / ((knots[i+n-1] - knots[i-1]) * (knots[i+n-1] - knots[i]));
			double const Q1 = - Q0 - Q2;

			return Q0 * coefs[i] + Q1 * coefs[i - 1] + Q2 * coefs[i - 2];
		}

		constexpr int nk = k > n ? k : k + 1;
		constexpr int nj = k > n ? j : j + 1;

		double const Akj = (x - knots[l - j]) / (knots[l - j + k] - knots[l - j]);
		double const Bkj = 1.0 - Akj;

		return Akj * der2_sum<nk,j>(l, x) + Bkj * der2_sum<nk,nj>(l, x);
	}

	inline double der2_sum(int l, double x) const
	{
		assert(n >= 2);
		return der2_sum<1,0>(l, x) * n * (n - 1);
	}

	/*
	 * for k = n-2:
	 * $S=\frac{d^{3}}{dx^{3}}N_{i}^{n}=Q_{0,i}N_{i}^{n-3}+Q_{1,i}N_{i+1}^{n-3}+Q_{2,i}N_{i+2}^{n-3}+Q_{3,i}N_{i+3}^{n-3}$
	 * with 
	 * $Q_{0,i}=\frac{n}{u_{i+n}-u_{i}}\cdot\frac{n-1}{u_{i+n-1}-u_{i}}\cdot\frac{n-2}{u_{i+n-2}-u_{i}}$
	 * $Q_{1,i-1}=-\frac{n}{u_{i+n}-u_{i}}\cdot\frac{n-1}{u_{i+n-1}-u_{i}}\cdot\frac{n-2}{u_{i+n-2}-u_{i}}-\frac{n}{u_{i+n-1}-u_{i-1}}\cdot\frac{n-1}{u_{i+n-2}-u_{i-1}}\cdot\frac{n-2}{u_{i+n-2}-u_{i}}-\frac{n}{u_{i+n-1}-u_{i-1}}\cdot\frac{n-1}{u_{i+n-1}-u_{i}}\cdot\frac{n-2}{u_{i+n-2}-u_{i}}$
	 * $Q_{2,i-2}=\frac{n}{u_{i+n-2}-u_{i-2}}\cdot\frac{n-1}{u_{i+n-2}-u_{i-1}}\cdot\frac{n-2}{u_{i+n-2}-u_{i}}+\frac{n}{u_{i+n-1}-u_{i-1}}\cdot\frac{n-1}{u_{i+n-2}-u_{i-1}}\cdot\frac{n-2}{u_{i+n-2}-u_{i}}+\frac{n}{u_{i+n-1}-u_{i-1}}\cdot\frac{n-1}{u_{i+n-1}-u_{i}}\cdot\frac{n-2}{u_{i+n-2}-u_{i}}$
	 * $Q_{3,i-3}=-\frac{n}{u_{i+n-2}-u_{i-2}}\cdot\frac{n-1}{u_{i+n-2}-u_{i-1}}\cdot\frac{n-2}{u_{i+n-2}-u_{i}}$
	 *
	 */
	template <int k, int j>
	inline double der3_sum(int l, double x) const
	{
		if (k == n - 2)
		{
			int const i = l - j;

			double const Q0 = 
				1.0 / ((knots[i+n] - knots[i]) * (knots[i+n-1] - knots[i]) * (knots[i+n-2] - knots[i]));
			double const Q3 = 
				-1.0 / ((knots[i+n-2] - knots[i-2]) * (knots[i+n-2] - knots[i-1]) * (knots[i+n-2] - knots[i]));
			double const Q1 = 
				-Q0
				-1.0 / ((knots[i+n-1] - knots[i-1]) * (knots[i+n-2] - knots[i-1]) * (knots[i+n-2] - knots[i]))
				-1.0 / ((knots[i+n-1] - knots[i-1]) * (knots[i+n-1] - knots[i]) * (knots[i+n-2] - knots[i]));
			double const Q2 = 
				-Q3
				+1.0 / ((knots[i+n-1] - knots[i-1]) * (knots[i+n-2] - knots[i-1]) * (knots[i+n-2] - knots[i]))
				+1.0 / ((knots[i+n-1] - knots[i-1]) * (knots[i+n-1] - knots[i]) * (knots[i+n-2] - knots[i]));

			return Q0 * coefs[i] + Q1 * coefs[i - 1] + Q2 * coefs[i - 2] + Q3 * coefs[i - 3];
		}

		constexpr int nk = k > n ? k : k + 1;
		constexpr int nj = k > n ? j : j + 1;

		double const Akj = (x - knots[l - j]) / (knots[l - j + k] - knots[l - j]);
		double const Bkj = 1.0 - Akj;

		return Akj * der3_sum<nk,j>(l, x) + Bkj * der3_sum<nk,nj>(l, x);
	}

	inline double der3_sum(int l, double x) const
	{
		assert(n >= 3);
		return der3_sum<1,0>(l, x) * n * (n - 1) * (n - 2);
	}

	inline int get_knot(double x) const
	{
		auto const& i = std::upper_bound(knots.begin(), knots.end(), x);
		return std::distance(knots.begin(), i - 1);
	}

public:
	non_uniform_spline(
		std::vector<double> const& knots, 
		std::vector<double> const& coefs 
		) : knots(knots), coefs(coefs)
	{
	}

	~non_uniform_spline()
	{
	}

	virtual double val(double x) const
	{
		int const l = get_knot(x);
		return sum(l, x);
	}

	virtual double der(double x) const
	{
		int const l = get_knot(x);
		return der_sum(l, x);
	}

	virtual double der2(double x) const
	{
		int const l = get_knot(x);
		return der2_sum(l, x);
	}

	virtual double der3(double x) const
	{
		int const l = get_knot(x);
		return der3_sum(l, x);
	}
};


/*
 * uniform knots spline
 */
template <int n /* degree of the spline */>
class uniform_spline : 
	public basic_spline
{
private:
	double 	min, max, step;
	std::vector<double> coefs;

private:
	template <int k, int j>
	inline double sum(int l, double x) const
	{
		if (k == n + 1)
			return coefs[l - j];

		constexpr int nk = k > n ? k : k + 1;
		constexpr int nj = k > n ? j : j + 1;

		double const Akj = (x - l + j);
		double const Bkj = -(x - l + j -k);

		return Akj * sum<nk,j>(l, x) + Bkj * sum<nk,nj>(l, x);
	}

	inline double sum(int l, double x) const
	{
		return sum<1,0>(l, (x - min) / step) / fact(n);
	}

	template <int k, int j>
	inline double der_sum(int l, double x) const
	{
		if (k == n)
			return coefs[l - j] - coefs[l - j - 1];

		constexpr int nk = k > n ? k : k + 1;
		constexpr int nj = k > n ? j : j + 1;

		double const Akj = (x - l + j);
		double const Bkj = -(x - l + j -k);

		return Akj * der_sum<nk,j>(l, x) + Bkj * der_sum<nk,nj>(l, x);
	}

	inline double der_sum(int l, double x) const
	{
		assert(n >= 1);
		return der_sum<1,0>(l, (x - min) / step) / (fact(n - 1) * step);
	}

	template <int k, int j>
	inline double der2_sum(int l, double x) const
	{
		if (k == n - 1)
			return coefs[l-j] - 2.0 * coefs[l-j-1] + coefs[l-j-2];

		constexpr int nk = k > n ? k : k + 1;
		constexpr int nj = k > n ? j : j + 1;

		double const Akj = (x - l + j);
		double const Bkj = -(x - l + j -k);

		return Akj * der2_sum<nk,j>(l, x) + Bkj * der2_sum<nk,nj>(l, x);
	}

	inline double der2_sum(int l, double x) const
	{
		assert(n >= 2);
		return der2_sum<1,0>(l, (x - min) / step) / (fact(n-2) * step * step);
	}

	template <int k, int j>
	inline double der3_sum(int l, double x) const
	{
		if (k == n - 2)
			return coefs[l-j] - 3.0 * coefs[l-j-1] + 3.0 * coefs[l-j-2] - coefs[l-j-3];

		constexpr int nk = k > n ? k : k + 1;
		constexpr int nj = k > n ? j : j + 1;

		double const Akj = (x - l + j);
		double const Bkj = -(x - l + j -k);

		return Akj * der3_sum<nk,j>(l, x) + Bkj * der3_sum<nk,nj>(l, x);
	}

	inline double der3_sum(int l, double x) const
	{
		assert(n >= 3);
		return der3_sum<1,0>(l, (x - min) / step) / (fact(n-3) * step * step * step);
	}

	inline int get_knot(double x) const
	{
		return static_cast<int>((x - min) / step);
	}

public:
	uniform_spline(double min, double max, std::vector<double> const& coefs) : 
		min(min), max(max), coefs(coefs)
	{
		step = (max - min) / (coefs.size() - 1);
	}

	~uniform_spline()
	{
	}

	virtual double val(double x) const
	{
		int const l = get_knot(x);
		return sum(l, x);
	}

	virtual double der(double x) const
	{
		int const l = get_knot(x);
		return der_sum(l, x);
	}

	virtual double der2(double x) const
	{
		int const l = get_knot(x);
		return der2_sum(l, x);
	}

	virtual double der3(double x) const
	{
		int const l = get_knot(x);
		return der3_sum(l, x);
	}

};


/*
 * spline wrapper
 */
bool is_uniform(std::vector<double> const& v, double eps)
{
	assert(v.size() > 1);

	double const max = *v.rbegin();
	double const min = *v.begin();
	double const step = (max - min) / (v.size() - 1);

	for (size_t i = 1; i < v.size(); ++ i)
	{
		if (fabs(v[i] - v[i-1] - step) >= eps)
			return false;
	}

	return true;
}

spline::spline(
	unsigned int degree, 
	std::vector<double> const& knots, 
	std::vector<double> const& coefs, 
	spline_extrapolation extrapolation
) : 
	extrapolation(extrapolation),
	min(knots[degree]),
	max(knots[knots.size() - degree - 1]),
	period(max - min)
{
	if (knots.size() <= degree * 2)
		throw std::invalid_argument("parameter knots is invalid");

	if (is_uniform(knots, 1e-10))
	{
		double const fisrst = *knots.begin();
		double const last = *knots.rbegin();

		switch (degree)
		{
		case 1: s.reset(new uniform_spline<1>(fisrst, last, coefs)); break;
		case 2: s.reset(new uniform_spline<2>(fisrst, last, coefs)); break;
		case 3: s.reset(new uniform_spline<3>(fisrst, last, coefs)); break;
		case 4: s.reset(new uniform_spline<4>(fisrst, last, coefs)); break;
		case 5: s.reset(new uniform_spline<5>(fisrst, last, coefs)); break;
		default: throw std::runtime_error("unsupported degree of spline: " + std::to_string(degree));
		}

		std::cout << "uniform spline created\n";
	}
	else
	{
		switch (degree)
		{
		case 1: s.reset(new non_uniform_spline<1>(knots, coefs)); break;
		case 2: s.reset(new non_uniform_spline<2>(knots, coefs)); break;
		case 3: s.reset(new non_uniform_spline<3>(knots, coefs)); break;
		case 4: s.reset(new non_uniform_spline<4>(knots, coefs)); break;
		case 5: s.reset(new non_uniform_spline<5>(knots, coefs)); break;
		default: throw std::runtime_error("unsupported degree of spline: " + std::to_string(degree));
		}

		std::cout << "non-uniform spline created\n";
	}
}

inline double spline::fix_arg(double x) const
{
	switch (extrapolation)
	{
	case spline_extrapolation::none:
		if (x < min || x > max) throw std::invalid_argument("the value is out of range");
		return x;
	case spline_extrapolation::periodic: return x - period * floor((x - min) / period);
	case spline_extrapolation::end_value: return std::max(min, std::min(max, x));
	default: throw std::invalid_argument("unknown extrapolation method");
	}
}

double spline::val(double x) const
{
	x = fix_arg(x);
	return s->val(x);
}

double spline::der(double x) const
{
	x = fix_arg(x);
	return s->der(x);
}

double spline::der2(double x) const
{
	x = fix_arg(x);
	return s->der2(x);
}

double spline::der3(double x) const
{
	x = fix_arg(x);
	return s->der3(x);
}
