splines
=======

the class evaluates value of spline function and its derivatives

```c++
spline::spline(
	unsigned int degree, // degree of the spline
	std::vector<double> const& knots, // array of knots
	std::vector<double> const& coefs, // array of coefficients
	spline_extrapolation extrapolation // type of extrapolation: none | periodic | end_value
)

spline::val(double x) // evaluates value of spline
spline::der(double x) // evaluates derivative of spline
spline::der2(double x) // evaluates 2d derivative of spline
spline::der3(double x) // evaluates 3d derivative of spline
```

The *coefficients* and *knots* parameters of the spline can be evaluated using *scipy.interpolate.spline* function of 
*python* or *spapi* of *MatLab*.
