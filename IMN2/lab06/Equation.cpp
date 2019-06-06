//
// Created by Kristina on 30-Oct-18.
//
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <zmouse.h>
#include "Equation.h"

Equation::Equation(int nx, int ny, double xp, double xk, double yp, double yk)
        :min(xp, yp, 1), max(xk, yk, nx*ny-1), N(nx*ny), E(N,N), O(N,N)
{
    dx = (max.x - min.x) / (nx-1);
    dy = (max.y - min.y) / (ny-1);
    // elements_size = 2*(nx-1)*(ny-1);

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
            E(i-1,j) = 0.; E(j,i-1) = 0.;
            O(i-1,j) = 0.; O(j,i-1) = 0.;
        }
        E(i-1,i-1) = 2000.;
        O(i-1,i-1) = 1.;
    }
}




void Equation::initialize() {
    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            E(i,j) = 0.;
            O(i,j) = 0.;
        }
    }
    for (auto &element : elements) {
        element.fillLocalStiffnessMatrix();
        element.fillLocalMassMatrix();

        element.fillGlobalStiffnessMatrix(E);
        element.fillGlobalMassMatrix(O);
    }
    setBoundaries();
//    std::cout << "\n\nMatrix E:\n" << E << "\n";
//    std::cout << "\n\nMatrix O:\n" << O << "\n";
//    printMatrix(E,26,26);
//    printMatrix(O,26,26);
}

void Equation::solve() {
    initialize();
    solver.compute(E,O);
    //printVector(solver.eigenvalues(), N);
    //std::cout << "\nWartosci wlasne:\n" << solver.eigenvalues() << "\n";

}


void Equation::computeU(std::ofstream &file, int k)  // k - indeks wartosci wlasnej, [0,N)
{
    double x , y;
    double dx0, dy0;
    int n = 100;
    dx0 = (max.x - min.x)/n; dy0 = (max.y - min.y)/n;
    x = min.x; y = min.y;

    c = solver.eigenvectors().col(k);
    std::cout << "\n\nE = " << solver.eigenvalues()[k]  << "\n";

    while(y <= max.y) {
        while (x <= max.x) {
            for (auto &element : elements) {
                //std::cout << m << "\n";
                if (element.isInside(x, y)) {
                    //std::cout << "element = " <<  m << "  , x = " << x << ", y = " << y << "\n";
                    double r = element.computeUValue(x, y, c);
                    u.push_back(r);
                    break;
                }
            }
            file << x << "\t" << y << "\t" << u[u.size()-1] << "\n";
            //std::cout << u[u.size() - 1] << "\n";
            x += dx0;
        }
        file << "\n";
        x = min.x;
        y += dy0;
    }
}

void Equation::computeU(std::ofstream &file, Eigen::VectorXd v)
{
    double x , y;
    double dx0, dy0;
    int n = 100;
    dx0 = (max.x - min.x)/n; dy0 = (max.y - min.y)/n;
    x = min.x; y = min.y;

    while(x <= max.x) {
        while (y <= max.y) {
            for (auto &element : elements) {
                //std::cout << m << "\n";
                if (element.isInside(x, y)) {
                    //std::cout << "element = " <<  m << "  , x = " << x << ", y = " << y << "\n";
                    double r = element.computeUValue(x, y, v);
                    u.push_back(r);
                    break;
                }
            }
            file << x << " " << y << " " << u[u.size()-1] << "\n";
            //std::cout << u[u.size() - 1] << "\n";
            y += dy0;
        }
        file << "\n";
        y = min.y;
        x += dx0;
    }
}



void Equation::printMatrix(Eigen::MatrixXd A, int n_rows, int n_cols) {
    std::cout << "Macierz :\n";
    for (int i = 0; i < n_rows ; ++i) {
        std::cout << "[" << i << "]\t";
        for (int j = 0; j < n_cols; ++j) {
            std::cout << A(i,j) << "\t";
        }
        std::cout << "\n";
    }
}

void Equation::printVector(Eigen::VectorXd v, int n_rows) {
    std::cout << "Wektor :\n";
    for (int i = 0; i < n_rows; ++i) {
        std::cout << "[" << i << "]\t" << v(i) << "\n";
    }
}

void Equation::writeMesh(std::ofstream &file) {
    for (auto &element : elements) {
        for (int l = 1; l <=3 ; ++l) {
            file << element.get(l).x << "\t" << element.get(l).y << "\t" << element.get(l).global_n << "\n";
        }
        file << element.get(1).x << "\t" << element.get(1).y << "\t" << element.get(1).global_n << "\n";
        file << "\n";
    }
}

