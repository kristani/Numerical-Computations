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

const double PI = acos(-1.);

class Equation {
public:
    Equation(double (*_ro)(double,double), int nx, int ny, double xp=0., double xk=PI, double yp=0., double yk=PI);

    void solve();
    double computeFunctionalIntegral();

    void computeU(std::ofstream& file);

    void writeMesh(std::ofstream& file);

    void printMatrix(int n_rows, int n_cols);
    void printVector(int n_rows);
private:


    int M, N;
    Node min, max;
    double dx, dy;
    std::vector<Element> elements;
    std::set<Node> boundaryNodes;
    std::vector<double> u;

    double (*ro)(double,double);
    Eigen::MatrixXd S;

    Eigen::VectorXd F, c;

    void findBoundaryNodes();
    void initialize();
    void setBoundary();
};


#endif //LAB04_EQUATION_H
