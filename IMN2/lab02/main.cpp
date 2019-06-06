#include "Equation.hpp"

double baza1(int i, double x)
{
	return cos( (i-1./2)*x )*exp(-x);
}

double baza2(int i, double x)
{
	return (x-PI)*(x+PI)*pow(x, i-1);
}

double f(double x)
{
	return 2*(1-3*x+x*x-PI*PI)*exp(-x);
}

double u_d(double x)
{
	return (x-PI)*(x+PI)*exp(-x);
}

int main()
{
	std::ofstream file;
	char numer[20];
	DiffOperator L(1., -1.);

	std::string nazwa = "zad1/kolokacji_";

	for(int N=6; N<=10; N++)
	{
		sprintf(numer, "%d.txt", N);
		file.open(nazwa+numer);

		Equation e(L, f, baza1, N, -PI, PI, u_d);
		e.solve(file,1);
	}


	nazwa = "zad2/kwadratow_";
	for(int N=6; N<=10; N++)
	{
		sprintf(numer, "%d.txt", N);
		file.open(nazwa+numer);

		Equation e(L, f, baza1, N, -PI, PI, u_d);
		e.solve(file,2);
	}

	nazwa = "zad3/galerkin_";
	for(int N=6; N<=10; N++)
	{
		sprintf(numer, "%d.txt", N);
		file.open(nazwa+numer);

		Equation e(L, f, baza1, N, -PI, PI, u_d);
		e.solve(file,3);
	}


	/************************** Zad4 **********************/
	nazwa = "zad4/kolokacji_";

	for(int N=6; N<=10; N++)
	{
		sprintf(numer, "%d.txt", N);
		file.open(nazwa+numer);

		Equation e(L, f, baza2, N, -PI, PI, u_d);
		e.solve(file,1);
	}


	nazwa = "zad4/kwadratow_";
	for(int N=6; N<=10; N++)
	{
		sprintf(numer, "%d.txt", N);
		file.open(nazwa+numer);

		Equation e(L, f, baza2, N, -PI, PI, u_d);
		e.solve(file,2);
	}

	nazwa = "zad4/galerkin_";
	for(int N=6; N<=10; N++)
	{
		sprintf(numer, "%d.txt", N);
		file.open(nazwa+numer);

		Equation e(L, f, baza2, N, -PI, PI, u_d);
		e.solve(file,3);
	}



	return 0;
}