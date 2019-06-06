//
// Created by Kristina on 30-Oct-18.
//

#ifndef LAB04_EQUATION_H
#define LAB04_EQUATION_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <set>

#include "../Eigen/Eigenvalues"
#include "../Eigen/Dense"

#include "Element.h"
#include "CircumCircle.h"

//const double PI = acos(-1.);

class Equation {
public:
    Equation(double (*_ro)(double,double), int nx, int ny, double xp, double xk, double yp, double yk);

    void check()
    {
        std::cout << "isInside: " << elements[0].isInside(-1.5, -1.9) << "\n";
    }

    void solve();
    void computeU(std::ofstream &file) {
        double x , y;
        double dx0, dy0;
        int n = 100;
        dx0 = (max.x - min.x)/n; dy0 = (max.y - min.y)/n;
        x = min.x; y = min.y;

        while(y <= max.y) {
            while (x <= max.x) {
                for (int m = 0; m < M; ++m) {
                    //std::cout << m << "\n";
                    if (elements[m].isInside(x, y)) {
                        //std::cout << "element = " <<  m << "  , x = " << x << ", y = " << y << "\n";
                        double r = elements[m].computeUValue(x, y, c);
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
    void printMatrix(int n_rows, int n_cols);
    void printVector(int n_rows);
    void writeMesh(std::ofstream &file)
    {
        for (auto &element : elements) {
            for (int l = 1; l <=3 ; ++l) {
                file << element.get(l).x << "\t" << element.get(l).y << "\t" << element.get(l).global_n << "\n";
            }
            file << element.get(1).x << "\t" << element.get(1).y << "\t" << element.get(1).global_n << "\n";
            file << "\n";
        }
    }
private:
    int M, N;

    Node min, max;
    double dx, dy;
    std::vector<Node> nodes;

    std::vector<Element> elements;
    std::set<Node> boundaryNodes;
    std::vector<double> u;
    double (*ro)(double,double);

    Eigen::MatrixXd S;
    Eigen::VectorXd F, c;

    void createMesh(int nx);
    std::vector<int> findInnerNodes();

    void initialize();
    void setBoundaries();

};


#endif //LAB04_EQUATION_H
