#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>

#include "../Eigen/Eigenvalues"
#include "../Eigen/Dense"

#include "../Gauss_Legendre.hpp"

const double PI = acos(-1);

// klasa reprezentujaca rownanie wlasne
class Equation
{
public:
	Equation()=default;
	Equation(double _xmax, int _M, double _alpha);
	~Equation();

	double basis(int, double) const;
	double diff_basis(int, double) const;
	double x_value(int m, double ksi) const;
	double ksi_value(int m, double x0) const;

	double integrate(int flag, int m, int i, int j);

	void initialize();
	void solve();
	void sort();

	void print(const Eigen::MatrixXd& A) const;

	void write_vector(std::ofstream& file, int mi) const;	
	void write_eigenvalues(std::ofstream& file) const;
	void write_u(std::ofstream& file, int) const;


private:

	double x_max;
	int M, N;
	double alpha;
	double dx;

	std::vector<double> x;
	std::vector<double> Jacobian;

	std::vector<double> E;
	std::vector<double> order;

	std::vector<double> u;

	Eigen::MatrixXd S, O, c;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;
};