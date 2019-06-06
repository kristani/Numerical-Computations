#include "Equation.hpp"


Equation::Equation(double _xmax, int _M, double _alpha):
					x_max(_xmax), M(_M), alpha(_alpha), N(2*M+1), E(N), order(N), S(N,N), O(N,N)
	{
		int sign;
		for(int k=1; k<=N; k++)
		{
			if( (2.*k-N-1.)/N < 0. ) sign=-1;
			else if( (2.*k-N-1.)/N == 0. ) sign = 0;
			else sign = 1;
			x.push_back( x_max * pow( abs(2.*k-N-1.)/(double)N ,alpha) * sign );
			//std::cout << x[k-1] << "\n";
		}
		//std::cout << N << " Jm \n";
		for(int m=1; m<=M; m++)
		{
			Jacobian.push_back( (x[2*m]-x[2*m-2])/2. );
			//std::cout << Jacobian[m-1] << "\n";
		}

		for(int i=0; i<N; i++)
		{
			for(int j=0; j<N; j++)
			{
				S(i,j) = 0.;
				O(i,j) = 0.;
			}
		}
	}


double  Equation::basis(int i, double ksi) const
{
	if(i==1) return ksi*(ksi-1.)/2.;
	else if(i==2) return -(ksi+1.)*(ksi-1.);
	else if(i==3) return ksi*(ksi+1.)/2.;
}


double Equation::diff_basis(int i, double ksi) const
{
	double d_ksi = 0.001;
	return ( basis(i,ksi+d_ksi) - basis(i,ksi-d_ksi) )/2./d_ksi;
}


double Equation::x_value(int m, double ksi) const
{
	return (x[2*m] + x[2*m-2])/2. + ksi*(x[2*m] - x[2*m-2])/2.;
}

double Equation::ksi_value(int m, double x0) const
{
	return 2*(x0 - (x[2*m] + x[2*m-2])/2.) / (x[2*m] - x[2*m-2]);
}


double Equation::integrate(int flag ,int m, int i, int j)
{
	const int n=5;
	Rosetta::GaussLegendreQuadrature<n> gl5;
	//gl5.print_roots_and_weights(std::cout);
	std::vector<double> weights = gl5.get_weights();
	std::vector<double> ksi = gl5.get_roots();
	double result = 0.;
	if(flag == 0){
		for(int k=0; k<n; k++)
		{
			result += weights[k] * (diff_basis(i,ksi[k]) * diff_basis(j,ksi[k]) / (2 * Jacobian[m-1]) + 
						Jacobian[m-1] * basis(i,ksi[k]) * basis(j,ksi[k]) * pow( x_value(m, ksi[k]) ,2)/2.);
		}
	}
	else if(flag == 1){
		for(int k=0; k<n; k++)
		{
			result += weights[k] * Jacobian[m-1] * basis(i,ksi[k]) * basis(j,ksi[k]);
		}
	}
	return result;
}


void Equation::initialize()
{
	int p, q;
	for(int m=1; m<=M; m++)
	{
		for(int i=1; i<=3; i++)
		{
			p = 2*m + (i-2);
			for(int j=1; j<=3; j++)
			{
				q = 2*m + (j-2);
				S(p-1,q-1) += integrate(0, m, i, j);
				O(p-1,q-1) += integrate(1, m, i, j);
			}
		}
	}
}



void Equation::solve()
{
	initialize();
	// std::cout << "\nmacierz S\n";
	// print(S);
	// std::cout << "\nmacierz O\n";
	// print(O);

	solver.compute(S,O);
	for(int i=0; i<N; i++) E[i] = solver.eigenvalues()[i]; 
	//std::cout << solver.eigenvalues();

	// std::cout << "\nEigen Vectors:\n";
	// std::cout << solver.eigenvectors().col(0);

}


void Equation::print(const Eigen::MatrixXd& A) const
{
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			std::cout << A(i,j) << "\t";
		}
		std::cout << "\n";
	}
}


void Equation::write_vector(std::ofstream& file, int mi) const
{
	for(int j=0; j<N; j++)
	{
		file << x[j] << "\t" << solver.eigenvectors().col(mi)[j] << "\n";
	}
	file.close();
}

void Equation::write_eigenvalues(std::ofstream& file) const
{
	file << alpha << "\t" << E[0] << "\t" << E[1] << "\t" << E[2] << "\t" << E[3] << "\t" << E[4] << "\n";
}

void Equation::write_u(std::ofstream& file, int mi) const
{
	int k=0, m=1, p;

	double x0 = x[0], ksi0, result;
	double dx = (x[N-1]-x[0])/500.;

	while( x0 <= x[N-1] )
	{
		if(x0 > x[k+2]){
			k += 2;
			m++;
		}
		ksi0 = ksi_value(m, x0);
		result = 0.;
		for(int i=1; i<=3; i++)
		{
			p = 2*m + (i-2);
			result += solver.eigenvectors().col(mi)[p-1] * basis(i,ksi0);
		}
		file << x0 << "\t" << result << "\n";
		x0 += dx;
	}
	file.close();
}

Equation::~Equation()
{
}