#include "Equation.hpp"
#include "nr/nr.h"
#include "nr/nrutil.h"
#include "nr/nrutil.c"
#include "nr/gaussj.c"
#include "nr/ludcmp.c"
#include "nr/lubksb.c"
#include "nr/gauleg.c"


Equation::Equation(DiffOperator _L, double (*_f)(double), double (*_basis)(int, double), int _N, double _xmin, double _xmax, double (*_u_dok)(double)):
					L{_L}, xmin(_xmin), xmax(_xmax), N(_N), basis(_basis), f(_f), u_dok(_u_dok)
	{
		dx = 2*PI/(N+1);
		K = N+2;
		for(int i=0; i<K; i++) x.push_back(xmin+i*dx);
		A = matrix(1,N,1,N), b = matrix(1,N,1,1);
		for(int i=1; i<=N; i++)
		{
			b[i][1] = 0.;
			for(int j=1; j<=N; j++)
			{
				A[i][j] = 0.;
			}
		}
	}

	void Equation::init_kolokacji()
	{
		for(int i=1; i<=N; i++)
		{
			b[i][1] = f(x[i]);
			for(int j=1; j<=N; j++)
			{
				A[i][j] = L.differentiate( basis,j, x[i] );
			}
		}
	}


	void Equation::init_kwadratow()
	{
		int m=40;
		float *part, *wagi;
		part = vector(1,m), wagi = vector(1,m);

		gauleg(xmin, xmax,part,wagi,m);

		for(int i=1; i<=N; i++)
		{
			for(int k=1; k<=m; k++) {b[i][1] += wagi[k] * f(part[k]) * L.differentiate( basis,i, (part[k]) );}
			for(int j=1; j<=N; j++)
			{
				for( int k=1; k<=m; k++)
				{
					A[i][j] += wagi[k] * L.differentiate( basis,i, (part[k]) ) * L.differentiate( basis,j, (part[k]) );
				}

			}
		}
		free_vector(part,1,m);
		free_vector(wagi,1,m);
	}

	void Equation::init_galerkin()
	{
		int m=40;
		float *part, *wagi;
		part = vector(1,m), wagi = vector(1,m);

		gauleg(xmin, xmax,part,wagi,m);

		for(int i=1; i<=N; i++)
		{
			for(int k=1; k<=m; k++) {b[i][1] += wagi[k] * f(part[k]) * basis(i, part[k]);}
			for(int j=1; j<=N; j++)
			{
				for( int k=1; k<=m; k++)
				{
					A[i][j] += wagi[k] * basis(i,part[k]) * L.differentiate( basis,j, (part[k]) );
				}

			}
		}
		free_vector(part,1,m);
		free_vector(wagi,1,m);
	}


void Equation::solve(std::ofstream& file,int i)
{
	//inicjalizacja
	if(i==1) init_kolokacji();
	else if(i==2) init_kwadratow();
	else if(i==3) init_galerkin();

	// std::cout << "init\n";
	// 	for(int j=1; j<=N; j++)
	// 	{
	// 		std::cout << b[j][1] << "\n";
	// 	}
	// 	std::cout << "\nA\n";
	// 	for(int i=1; i<=N; i++)
	// {
	// 	for(int j=1; j<=N; j++)
	// 	{
	// 		std::cout << A[i][j] << "\t";
	// 	}
	// 	std::cout << "\n";
	// }

	// rozwiazanie UR
	// std::cout << "\n\nhal\n";
	gaussj(A, N, b, 1);
	// int *indx;
	// float *d;
	// ludcmp(A, N, indx, d);
	// lubksb(A, N, indx, b);

	// for(int j=1; j<=N; j++)
	// 	{
	// 		std::cout << b[j][1] << "\n";
	// 	}


	// obliczenie u
	u.clear();
	u.push_back(0.);
	for(int i=1; i<=n; i++)
	{
		u.push_back(0.);
		for(int j=1; j<=N; j++)
		{
			u[i] += b[j][1]*basis(j, xmin+dx2*i);
		}
	}
	u.push_back(0);

	// wypisanie
	write_u(file);

}

void Equation::write_u(std::ofstream& file)
{
	double x;
	for(int i=0; i<u.size(); i++)
	{
		x = xmin+dx2*i;
		file << x << "\t" << u_dok(x) - u[i] << "\n";
	}
	file.close();
}

Equation::~Equation()
{
	free_matrix(b,1,N,1,1);
	free_matrix(A,1,N,1,N);
}