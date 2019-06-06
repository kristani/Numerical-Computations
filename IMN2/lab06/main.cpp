#include <iostream>
#include "Equation.h"


int main() {
    std::ofstream file;

    int n = 10;
    double L = 5.;
    Equation e(n, n, -L, L, -L, L);
    file.open("../results/mesh.txt");
    e.writeMesh(file);
    file.close();

    e.solve();

    std::string name= "../results/zad1/u_";

    for (int k = 0; k < 10; ++k) {
        file.open(name + std::to_string(k) + ".txt");
        e.computeU(file, k);
        file.close();
    }

    e.newmark();

    //e.check();

    return 0;
}