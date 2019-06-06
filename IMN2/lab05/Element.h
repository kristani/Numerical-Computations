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
    Element(Node* n1, Node* n2, Node* n3, Node* _min, Node* _max);

    Node& get(int l);
    double computeArea(Node n1, Node n2, Node n3);
    bool isInside(double x, double y);
    void checkLocalIndexes();
    void computeParameters()
    {
        jacobian = Jacobian();
        area = computeArea(*nodes[0],*nodes[1],*nodes[2]);
    }
    void fillLocalMatrix();
    void fillLocalVector(double (*ro)(double, double));

    void fillGlobalMatrix(Eigen::MatrixXd& S);
    void fillGlobalVector(Eigen::VectorXd& F);

    double computeUValue(double x, double y, Eigen::VectorXd &c) {
        double result = 0.;
        int p;
        for(int l=0; l<3; l++)
        {
            p = nodes[l]->global_n - 1;
            double fi1 = fi(l, ksi1(x,y), ksi2(x,y));
            double cp = c(p);
            result += c(p) * fi(l, ksi1(x,y), ksi2(x,y));
        }
        return result;
    }
    void write(std::ofstream& file){}
    void print();

private:

    Node *min_global, *max_global;

    double area;
    std::vector<Node*> nodes;
    Eigen::Matrix3d S_local;
    Eigen::Vector3d F_local;

    double jacobian;
    const static int n=7;


    double fi(int i, double ksi1, double ksi2);
    double d_fi(int flag, int i, double ksi1, double ksi2);
    double x(double ksi1, double ksi2);
    double y(double ksi1, double ksi2);

    double ksi1(double x, double y);
    double ksi2(double x, double y);

    double d_x(int flag, double ksi1, double ksi2);
    double d_y(int flag, double ksi1, double ksi2);
    double d_ksi1(int flag, double ksi1, double ksi2);
    double d_ksi2(int flag, double ksi1, double ksi2);
    double nabla_fi(int i, int j, double ksi1, double ksi2);
    double Jacobian();

    double integrate(int i, int j);   // for matrix
    double integrate(int j, double (*ro)(double, double));   // for vector

};


#endif //LAB04_ELEMENT_H
