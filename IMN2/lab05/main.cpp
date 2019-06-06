#include <iostream>
#include "Equation.h"

double ro(double x, double y)
{
    return exp(-(x*x + y*y)/2.);
}

int main() {
    std::ofstream file;

    int n = 10;
    Equation e(ro,n,n,-5.,5.,-5.,5.);
    file.open("../results/mesh.txt");
    e.writeMesh(file);
    file.close();

    e.solve();

    file.open("../results/u_10.txt");
    e.computeU(file);
    file.close();
    //e.check();

    return 0;
}