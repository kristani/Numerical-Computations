#include "Psi.hpp"
#include <iostream>

using namespace std;

int main()
{
	ofstream file;
	file.open("zad1.txt");
	Psi zad1;
	zad1.shooting(file);
	file.close();


	file.open("zad2.txt");
	Psi zad2;

	std::vector<double> zera = zad2.find_zeros();
	std::vector<double> zera_analit;
	zera_analit.push_back(pow(2.4048, 2)/2);
	zera_analit.push_back(pow(5.5200, 2)/2);
	zera_analit.push_back(pow(8.6537, 2)/2);
	zera_analit.push_back(pow(11.7915, 2)/2);

	file << "analityczne" << "\t" << "numeryczne" << "\n";
	for(int i=0; i<4; i++)
	{
		file << zera_analit[i] << "\t" << zera[i] << "\n";
	}	
	file.close();


	Psi zad3;
	zad3.compare("zad3");

	Psi zad4;
	zad3.compare("zad4", 1);

	return 0;

}