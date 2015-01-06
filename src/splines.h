#ifndef __splines_h__
#define __splines_h__

#include <vector>
#include <memory>


class basic_spline
{
public:
	virtual double val(double x) const = 0;
	virtual double der(double x) const = 0;
	virtual double der2(double x) const = 0;
	virtual double der3(double x) const = 0;
};

enum spline_extrapolation
{
	none,
	periodic,
	end_value
};

class spline
{
private:
	std::unique_ptr<basic_spline> 	s;
	spline_extrapolation const 		extrapolation;
	double const 					min, max, period;

	double fix_arg(double x) const;

public:
	spline(
		unsigned int degree, 
		std::vector<double> const& knots, 
		std::vector<double> const& coefs,
		spline_extrapolation extrapolation = spline_extrapolation::none
	);
	double val(double x) const;
	double der(double x) const;
	double der2(double x) const;
	double der3(double x) const;
};

#endif // __splines_h__
