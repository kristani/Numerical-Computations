#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <cstdio>


#include "nr.h"
#include "bessj0.c"
#include "bessj1.c"


class Psi
{
public:
	Psi(): R(N), U(N), g(N){}
	void initial(double r0, double r1)
	{
		R[0] = r0;
		R[1] = r1;
		U[0] = 0, U[1] = R[1]*sqrt(dr);
	}

	double compute_R(double E, double l=0)
	{
		double ri;
		for(int i=1; i<N-1; i++)
		{
			ri = dr * i;
			R[i+1] = ( R[i]*(2/dr/dr + l*l/(ri*ri) - 2*E) + R[i-1]*(-1/dr/dr + 1/(2*dr*ri)) ) / (1/dr/dr + 1/(2*dr*ri));
		}

		return R[N-1];
	}

	void shooting(std::ofstream& file, int l=0)
	{
		initial(1,1);
		double E = dE, R;
		
		while(E<=150)
		{
			R = compute_R(E);
			//std::cout << "jestem" <<E << "\t" << R << "\n";
			file << E << "\t" << R << "\n";

			E += dE;
		}
	}

	std::vector<double> find_zeros(int l=0)
	{
		if(l==0){initial(1,1);}
		else {initial(0,1);}

		double E = dE, R0 = 0., R1 = 0.;
		double E0 = 0., E1 =0., E_n;

		std::vector<double> zera;
		
		while(E<=150)
		{
			R1 = compute(E, l);
			
			//std::cout << "jestem" <<E << "\t" << R1 << "\n";
			
			// std::cout << "jestem " << E << "\t" << E-dE << "\t" << R0 << "\t" << R1 << "\n";
			if(R0*R1 < 0)
			{
				//std::cout << "jestem " << E << "\t" << E-dE << "\t" << R0 << "\t" << R1 << "\n";
				E0 = E-dE, E1 = E;
				//std::cout << "jestem " << E0 << "\t" << E1 << "\n";
				while( fabs(E1-E0) > 1e-6 )
				{
					E_n = E1 - R1*(E1 - E0)/(R1-R0);

					E0 = E1;
					E1 = E_n;
					R0 = R1, R1 = compute(E1, l);
				}
				//std::cout << "jestem tutaj " << E0 << "\t" << E1 << "\n";
				zera.push_back(E1);
			}


			R0 = compute(E, l);
			E += dE;
		}

		return move(zera);
	}

	void compare(std::string s0, int l=0)
	{
		std::vector<double> zera = find_zeros(l);
		std::vector<double> zera_analit = {2.4048, 5.5200, 8.6537, 11.7915};
		//std::vector<double> R_analit{N};
		for(int i=0; i<4; i++)
		{
			char name[20];
			sprintf(name, "_E%d.txt", i);
			std::string s(name);
			std::ofstream file(s0+s);

			compute(zera[i], l);
			for(int j=0; j<N; j++)
			{
				file << j*dr << "\t" << ( bessj0(zera_analit[i]*dr*j)-R[j] ) << "\n";
				//file << j*dr << "\t" << fabs(R[j]-R_analit[j]) << "\n";
			}
			std::cout << "jestem\n";
		}
		
	}

	double numerov(double E, double l)
	{
		for(int i=1; i<N; i++)
		{
			g[i] = 2*E + (1-4*l*l)/(4*dr*i*dr*i);
		}
		g[0] = g[1];
		for(int i=1; i<N-1; i++)
		{
			U[i+1] = ( U[i]*2*(1 - 5*dr*dr*g[i]/12) - U[i-1]*(1 + dr*dr*g[i-1]/12) ) / (1 + dr*dr*g[i+1]/12);
			R[i+1] = U[i+1]/sqrt(dr*i);
		}

		return R[N-1];
	}

	double compute(double E, double l=0)
	{
		if(l==0) {
			return compute_R(E);
		}
		else if(l==1)
		{
			return numerov(E,l);
		}
	}

private:
	const double L = 1.;
	const double dr = 0.01;
	const int N = 101; 
	
	std::vector<double> R;

	std::vector<double> U;
	std::vector<double> g;

	const double dE = 0.2; 
	const double E_k = 150;
	int l;		// moment pędu
};

