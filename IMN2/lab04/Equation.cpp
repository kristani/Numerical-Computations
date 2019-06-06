//
// Created by Kristina on 30-Oct-18.
//

#include "Equation.h"

Equation::Equation(double (*_ro)(double,double), int nx, int ny, double xp, double xk, double yp, double yk)
        : ro(_ro), min(xp, yp, 0), max(xk, yk, nx*ny), N(4*nx*ny), S(N,N), F(N)
{
    dx = (max.x - min.x) / (nx-1);
    dy = (max.y - min.y) / (ny-1);
    M = (nx-1)*(ny-1);

    for(int j = 1; j < nx; ++j)
    {
        for (int i = 1; i < ny; ++i)
        {
            elements.emplace_back(i, j, nx, dx, dy, min, max);
        }
    }

    for(int i=0; i<N; i++)
    {
        F(i) = 0.;
        for(int j=0; j<N; j++)
        {
            S(i,j) = 0.;
        }
    }
}

void Equation::initialize() {
    for(int m=1; m<=M; m++)
    {
        elements[m-1].addToGlobalMatrix(S);
        elements[m-1].addToGlobalVector(F, ro);
    }

    setBoundary();
}

void Equation::findBoundaryNodes() {
    std::vector<Node> temp;
    for (int m = 0; m < M; ++m) {
        temp = elements[m].getEdgeNodes();
        boundaryNodes.insert(temp.begin(), temp.end());
        temp.clear();
    }
}

void Equation::setBoundary() {
    findBoundaryNodes();

    int p ;
    for (const auto &node : boundaryNodes) {
        for (int i = 0; i < 3; ++i)
        {
            p = 4*(node.global_n - 1 ) + i;
            for (int j = 0; j < N; ++j) {S(p,j) = 0.; S(j,p) = 0.;}
            S(p,p) = 1.;
            F(p) = 0.;
        }
    }

}

void Equation::solve() {
    initialize();
    //printMatrix(N/2, N/4);
    //printVector(N/4);

    c = S.fullPivLu().solve(F);
}

double Equation::computeFunctionalIntegral() {
    double a =0.;
    for (int i = 0; i < N; ++i)
    {
        a += c(i)*F(i);
        for (int j = 0; j < N; ++j)
        {
            a -= c(i)*c(j)*S(i,j)/2.;
        }
    }
    //std::cout << a << '\n';
    return a;
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

void Equation::computeU(std::ofstream &file) {
    double x , y;
    double dx, dy;
    int n = 100;
    dx = (max.x - min.x)/n; dy = (max.y - min.y)/n;
    x = min.x; y = min.y;

    while(y <= max.y) {
        while (x <= max.x) {
            for (int m = 0; m < M; ++m) {
                //std::cout << m << "\n";
                if (elements[m].isElement(x, y)) {
                    //std::cout << "element = " <<  m << "  , x = " << x << ", y = " << y << "\n";
                    u.push_back(elements[m].computeUValue(x, y, c));
                    break;
                }
            }
            file << x << "\t" << y << "\t" << u[u.size() - 1] << "\n";
            //std::cout << u[u.size() - 1] << "\n";
            x += dx;
        }
        file << "\n";
        x = min.x;
        y += dy;
    }
}

void Equation::writeMesh(std::ofstream &file) {
    for (int m = 0; m < M; ++m)
    {
        for (int l = 1; l <=4 ; ++l) {
            file << elements[m].get(l).x << "\t" << elements[m].get(l).y << "\n";
        }
    }
}
