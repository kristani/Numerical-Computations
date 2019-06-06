#include "Equation.hpp"


double u_d(double x)
{
	return (x-PI)*(x+PI)*exp(-x);
}

int main()
{
	std::ofstream file;
	char numer[20];

	double x_m = 6.;
	int M;
	double alpha, d_alpha = 0.05;
	std::string nazwa;


	/**************** Wypisanie E(alpha) ********/
	M = 5;
	alpha = 0.4;
	file.open("results/M_5/E.txt");
	while(alpha <= 2.0 + d_alpha/10.)
	{
		Equation e(x_m, M, alpha);
		e.solve();
		e.write_eigenvalues(file);
		alpha += d_alpha;
	}
	file.close();


	M = 10;
	alpha = 0.4;
	file.open("results/M_10/E.txt");
	while(alpha <= 2.0 + d_alpha/10.)
	{
		Equation e(x_m, M, alpha);
		e.solve();
		e.write_eigenvalues(file);
		alpha += d_alpha;
	}
	file.close();


	M = 30;
	alpha = 0.4;
	file.open("results/M_30/E.txt");
	while(alpha <= 2.0 + d_alpha/10.)	
	{
		Equation e(x_m, M, alpha);
		e.solve();
		e.write_eigenvalues(file);
		alpha += d_alpha;
	}
	file.close();


	/************ Wypisywanie u(x) dla różnych M ***************/
	alpha = 1.4;
	M = 5;
	Equation e5(x_m, M, alpha);
	e5.solve();
	nazwa = "results/M_5/u_";
	for(int i=0; i<5; i++) 
	{
		sprintf(numer, "%d.txt", i);
		file.open(nazwa+numer);
		e5.write_u(file, i);
	}


	M = 10;
	Equation e10(x_m, M, alpha);
	e10.solve();
	nazwa = "results/M_10/u_";
	for(int i=0; i<5; i++) 
	{
		sprintf(numer, "%d.txt", i);
		file.open(nazwa+numer);
		e10.write_u(file, i);
	}


	M = 30;
	Equation e30(x_m, M, alpha);
	e30.solve();
	nazwa = "results/M_30/u_";
	for(int i=0; i<5; i++) 
	{
		sprintf(numer, "%d.txt", i);
		file.open(nazwa+numer);
		e30.write_u(file, i);
	}


	return 0;
}