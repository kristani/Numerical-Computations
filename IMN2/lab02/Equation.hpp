#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>

const double PI = acos(-1);

class DiffOperator
{
public:
	DiffOperator(){}
	DiffOperator(DiffOperator& L) = default;
	DiffOperator(double _a, double _b): a{_a}, b{_b}{}
	double differentiate(double (*v)(int, double), int i, double x)
	{
		double d1 = a*( v(i,x+dx) - 2*v(i,x) + v(i,x-dx) )/dx/dx;
		double d2 = b*( v(i,x+dx) - v(i,x-dx) )/dx/2;
		return d1+d2;
	}


private:
	double a;
	double b;
	const double dx = 0.001;

};

class Equation
{
public:
	Equation()=default;
	Equation(DiffOperator _L, double (*_f)(double), double (*_basis)(int, double), int _N, double _xmin, double _xmax, double (*_u_dok)(double));
	~Equation();

	void init_kolokacji();
	void init_kwadratow();
	void init_galerkin();

	void solve(std::ofstream& file, int i);

	void write_u(std::ofstream& file);

private:
	DiffOperator L;
	double xmin, xmax;
	int N, K;
	double dx;

	const int n = 200;
	const double dx2 = 2*PI/(n+1);

	std::vector<double> x;
	std::vector<double> u;
	double (*basis)(int, double);
	double (*f)(double);
	double (*u_dok)(double);

	float **A, **b, *c;
	
};