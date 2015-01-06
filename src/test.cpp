#include <iostream>
#include <fstream>
#include <memory.h>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <array>
#include "splines.h"


int ticks()
{
    timeval tv;
    gettimeofday(&tv, nullptr);
    return tv.tv_usec;
}

void read_array(std::istream& s, std::vector<double>& res)
{
	char c;

	s >> std::skipws >> c;
	if (c != '[')
		throw std::runtime_error("exptected '['");

	res.clear();

	while (true)
	{
		double x;

		s >> std::skipws >> x;
		res.push_back(x);
		s >> std::skipws >> c;

		if (c == ']')
			break;
		if (c != ',')
			throw std::runtime_error("exptected ',' or ']");
	}
}

struct spline_data
{
	int degree;
	std::vector<double> knots;
	std::vector<double> cp;
};

void read_spline_coefs(std::istream& s, spline_data& res)
{
	s >> std::skipws >> res.degree;
	read_array(s, res.knots);
	read_array(s, res.cp);
}

void read_spline_coefs(std::string const& filename, spline_data& res)
{
	std::ifstream f(filename, std::ifstream::in);
	if (!f.good()) throw std::runtime_error("can't open file");
	read_spline_coefs(f, res);	
}

int main()
{
	// spline
	spline_data data;
	read_spline_coefs("spline", data);
	spline s(data.degree, data.knots, data.cp, spline_extrapolation::periodic);

	std::vector<double> args;
	std::vector<double> res;

	// read
	while (!std::cin.eof())
	{
		double a;
		std::cin >> a;
		args.push_back(a);
	}

	// for (size_t i = 0; i < args.size(); ++ i)
	// 	std::cout << "sin(" << args[i] << ") = " << s.val(args[i]) << std::endl;

	// return 0;

	auto f = [](double x) { return sin(x + 1.2345); };
	auto df = [](double x) { return cos(x + 1.2345); };
	auto ddf = [](double x) { return -sin(x + 1.2345); };
	auto dddf = [](double x) { return -cos(x + 1.2345); };

	int t;
	int max_err_idx;

	/*
	 * evaluate function value
	 */ 
	res.clear();
	res.reserve(args.size());

	t = ticks();
	for (auto a : args)
		res.push_back(s.val(a));
	t = ticks() - t;

	max_err_idx = 0;
	for (size_t i = 0; i < args.size(); ++ i)
		max_err_idx = fabs(f(args[max_err_idx]) - res[max_err_idx]) < fabs(f(args[i]) - res[i]) ? 
			i : max_err_idx;

	std::cout << "max error: " << f(args[max_err_idx]) - res[max_err_idx] << " at " << args[max_err_idx] << std::endl;
	std::cout << "perf: " << t << "usec" << " for " << args.size() << " times" << std::endl;

	/*
	 * evaluate function derivative
	 */ 
	res.clear();
	res.reserve(args.size());

	t = ticks();
	for (auto a : args)
		res.push_back(s.der(a));
	t = ticks() - t;

	max_err_idx = 0;
	for (size_t i = 0; i < args.size(); ++ i)
		max_err_idx = fabs(df(args[max_err_idx]) - res[max_err_idx]) < fabs(df(args[i]) - res[i]) ? 
			i : max_err_idx;

	std::cout << "max error: " << df(args[max_err_idx]) - res[max_err_idx] << " at " << args[max_err_idx] << std::endl;
	std::cout << "perf: " << t << "usec" << " for " << args.size() << " times" << std::endl;

	/*
	 * evaluate function second derivative
	 */ 
	res.clear();
	res.reserve(args.size());

	t = ticks();
	for (auto a : args)
		res.push_back(s.der2(a));
	t = ticks() - t;

	max_err_idx = 0;
	for (size_t i = 0; i < args.size(); ++ i)
		max_err_idx = fabs(ddf(args[max_err_idx]) - res[max_err_idx]) < fabs(ddf(args[i]) - res[i]) ? 
			i : max_err_idx;

	std::cout << "max error: " << ddf(args[max_err_idx]) - res[max_err_idx] << " at " << args[max_err_idx] << std::endl;
	std::cout << "perf: " << t << "usec" << " for " << args.size() << " times" << std::endl;

	/*
	 * evaluate function third derivative
	 */ 
	res.clear();
	res.reserve(args.size());

	t = ticks();
	for (auto a : args)
		res.push_back(s.der3(a));
	t = ticks() - t;

	max_err_idx = 0;
	for (size_t i = 0; i < args.size(); ++ i)
		max_err_idx = fabs(dddf(args[max_err_idx]) - res[max_err_idx]) < fabs(dddf(args[i]) - res[i]) ? 
			i : max_err_idx;

	std::cout << "max error: " << dddf(args[max_err_idx]) - res[max_err_idx] << " at " << args[max_err_idx] << std::endl;
	std::cout << "perf: " << t << "usec" << " for " << args.size() << " times" << std::endl;

    return 0;
}

