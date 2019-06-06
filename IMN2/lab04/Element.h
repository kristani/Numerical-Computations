//
// Created by Kristina on 30-Oct-18.
//

#ifndef LAB04_ELEMENT_H
#define LAB04_ELEMENT_H

#include <vector>
#include <math.h>
#include <cmath>

#include "../Eigen/Eigenvalues"
#include "../Eigen/Dense"
#include "../Gauss_Legendre.hpp"

#include "Node.h"

class Element {
public:
    Element(int i, int j, int nx, double dx, double dy, Node _min, Node _max);

    Node& get(int l);
    void addToGlobalMatrix(Eigen::MatrixXd& S);
    void addToGlobalVector(Eigen::VectorXd& F, double (*ro)(double,double));
    std::vector<Node> getEdgeNodes();

    bool isElement(double x, double y);

    double ksi1(double x);
    double ksi2(double y);
    double computeUValue(double x, double y, Eigen::VectorXd& c);

private:

    Node min_global, max_global;
    Node min, max;
    std::vector<Node> nodes;

    const static int n=21;
    Rosetta::GaussLegendreQuadrature<n> gl5;

    int get_alpha(int l);
    int get_beta(int l);

    double fi(int a, int b, double ksi);        // Hermite basis functions
    double d_fi(int a, int b, double ksi);
    double Jacobian();
    double w(int i, double ksi1, double ksi2);
    double x(double ksi1, double ksi2);
    double y(double ksi1, double ksi2);

    double integrate(int l1, int i1, int i2, double (*ro)(double,double));   // for vector
    double integrate(int l1, int l2, int i1, int i2, int j1, int j2);   // for matrix
    void matrixForNodes(int l1, int l2, Eigen::MatrixXd& S_node);
};


#endif //LAB04_ELEMENT_H
