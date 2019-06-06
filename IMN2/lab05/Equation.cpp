//
// Created by Kristina on 30-Oct-18.
//
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "Equation.h"

Equation::Equation(double (*_ro)(double,double), int nx, int ny, double xp, double xk, double yp, double yk)
        : ro(_ro), min(xp, yp, 1), max(xk, yk, nx*ny-1), N(nx*ny), S(N,N), F(N)
{
    dx = (max.x - min.x) / (nx-1);
    dy = (max.y - min.y) / (ny-1);
    M = 2*(nx-1)*(ny-1);

    for(int i=0; i<N; i++)
    {
        F(i) = 0.;
        for(int j=0; j<N; j++)
        {
            S(i,j) = 0.;
        }
    }
    // tworzenie siatki
    createMesh(nx);

}

void Equation::createMesh(int nx) {
    // generacja wezlow
    for (int j = 0; j < nx; ++j) {
        for (int i = 0; i < nx; ++i) {
            nodes.emplace_back(min.x + i*dx, min.y + j*dy, nx*j+i+1);
        }
    }
    //    std::cout << dx << "\n";
    //    std::cout << nodes[4].x << " " << nodes[4].y << " " << nodes[4].global_n<< "\n";

    //losowe przesunięcie
    std::vector<int> inner = findInnerNodes();

    srand ( static_cast<unsigned int>(time(nullptr)) );
    std::vector<Node> displacements;
    auto rand_d = [](){return (rand()/(double)RAND_MAX - 0.5)/500.;};
    for (int k = 0; k < inner.size(); ++k) {
        displacements.emplace_back( rand_d(), rand_d(), 0 );
    }
    for(int k=0; k < inner.size(); k++)
    {
        nodes[ inner[k] ] += displacements[k];
    }
    //for(auto x:displacements) std::cout << x.x << "," << x.y << "\n";

    //triangulacja
    int n = nodes.size();
    for (int i1 = 0; i1 < n-2; ++i1) {
        int i2 = i1+1;
        while(i2 < n - 1)
        {
            int i3 = i2+1;
            while(i3 < n)
            {
                double w = (nodes[i1].x - nodes[i2].x)*(nodes[i1].y - nodes[i3].y) -
                        (nodes[i1].x - nodes[i3].x)*(nodes[i1].y - nodes[i2].y);
                if( w != 0 )
                {
                    CircumCircle circle(nodes[i1], nodes[i2], nodes[i3]);
                    bool flag = true;
                    for (int k = 0; k < n; ++k) {
                        if(k != i1 && k != i2 && k != i3) {
                            if(circle.isInside(nodes[k])) {
                                flag = false; break;
                            }
                        }
                    }
                    if(flag) elements.emplace_back(&nodes[i1], &nodes[i2], &nodes[i3], &nodes[0], &nodes[N-1]);

                }

                i3++;
            }
            i2++;
        }

    }

    for(int k=0; k < inner.size(); k++)
    {
        nodes[ inner[k] ] -= displacements[k];
    }
    std::cout << "koniec: elementow = " << elements.size() << "; wezlow = " << nodes.size() << "\n";
    for(auto &x:elements)
    {
        x.checkLocalIndexes();
        x.computeParameters();
    }
    //for(auto x:elements) x.print();
}

std::vector<int> Equation::findInnerNodes() {
    std::vector<int> inner;
    for (int i = 0; i < nodes.size(); ++i) {
        if( nodes[i].x != min.x && nodes[i].x != max.x && nodes[i].y != min.y && nodes[i].y != max.y )
            inner.push_back(i);
    }
    return std::move(inner);
}

void Equation::setBoundaries() {
    std::vector<int> boundaryNodesIndexes;
    for(auto& node : nodes)
    {
        if(node.x == nodes[0].x || node.y == nodes[0].y || node.x == nodes[nodes.size()-1].x || node.y == nodes[nodes.size()-1].y){
            boundaryNodesIndexes.push_back(node.global_n);
        }
    }

    for (auto i : boundaryNodesIndexes) {
        for (int j = 0; j < N; ++j) {
            S(i-1,j) = 0.;
            //S(j,i-1) = 0.;
        }
        S(i-1,i-1) = 1.;
        F(i-1) = 0.;
    }
}

void Equation::initialize() {
    for (auto &element : elements) {
        element.fillLocalMatrix();
        element.fillLocalVector(ro);

        element.fillGlobalMatrix(S);
        element.fillGlobalVector(F);
    }
    setBoundaries();
    //printMatrix(16,16);
    //printVector(16);
}

void Equation::solve() {
    initialize();
    printMatrix(24,24);
    printVector(24);
    c = S.fullPivLu().solve(F);
    std::cout << "\nWektor c:\n" << c << "\n";
}

void Equation::printMatrix(int n_rows, int n_cols) {
    std::cout << "Macierz S\n";
    for (int i = 0; i < n_rows ; ++i) {
        std::cout << "[" << i << "]\t";
        for (int j = 0; j < n_cols; ++j) {
            std::cout << S(i,j) << "\t";
        }
        std::cout << "\n";
    }
}

void Equation::printVector(int n_rows) {
    std::cout << "Wektor F:\n";
    for (int i = 0; i < n_rows; ++i) {
        std::cout << "[" << i << "]\t" << F(i) << "\n";
    }
}

