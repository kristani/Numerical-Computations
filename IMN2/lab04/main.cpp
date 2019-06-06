#include <iostream>
#include "Equation.h"

double ro(double x, double y)
{
    return sin(2*y)*sin(x)*sin(x);
}

int main() {
    int n;
    std::ofstream file;

    /************* Rysowanie siatki **************/
    n=3;
    Equation e(ro, n, n);
    file.open("../results/siatka.txt");
    e.writeMesh(file);
    file.close();
    /**********************************************/

    /************* Obliczanie całki funkcjonalnej **************/
    double a_dok = -0.1805396133;
    file.open("../results/a_integral.txt");
    file << "nx" << "\t" << "a_num" << "\t" << "a_dok" << "\n";

    e.solve();
    file << n << "\t" << e.computeFunctionalIntegral() << "\t" << a_dok << "\n";

    for(n=5; n<=20; n+=5)
    {
        Equation e0(ro, n, n);
        e0.solve();
        file << n << "\t" << e0.computeFunctionalIntegral() << "\t" << a_dok  << "\n";
        std::cout << n << " done\n";
    }
    file.close();
    /**************************************************************/

    /************* Rysowanie U **************/
    file.open("../results/u_nx_3.txt");
    e.computeU(file);
    file.close();

    n=10;
    Equation e2(ro, n, n);
    e2.solve();
    file.open("../results/u_nx_10.txt");
    e2.computeU(file);
    file.close();
    return 0;
}