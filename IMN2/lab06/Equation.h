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
    Equation(int nx, int ny, double xp, double xk, double yp, double yk);

    void check()
    {
        std::cout << "isInside: " << elements[0].isInside(-1.5, -1.9) << "\n";
    }

    void solve();
    void newmark()
    {
        Eigen::VectorXd c2 = solver.eigenvectors().col(1);
        Eigen::VectorXd c3 = solver.eigenvectors().col(2);
        double w2 = sqrt( solver.eigenvalues()[1] );
        double w3 = sqrt( solver.eigenvalues()[2] );

        Eigen::VectorXd y0 = c2;
        Eigen::VectorXd v0 = w3*c3;
        Eigen::VectorXd b, y, v;

        double beta = 1.4;
        double T = 2*M_PI/w2, dt = T/1e4, t = 0.;
        int it = 0;

        Eigen::MatrixXd A = O + beta*dt*dt*E;
        A_c.compute(A);
        O_c.compute(O);

        std::ofstream file1, file2, file3, file4, file_a;
        file1.open("../results/zad2/val_1.txt");
        file2.open("../results/zad2/val_2.txt");
        file3.open("../results/zad2/val_3.txt");
        file4.open("../results/zad2/val_4.txt");

        std::string name = "../results/zad3/u_";
        while(t <= T + dt/10.)
        {
            b = O*y0 + dt*O*v0 + (2.*beta - 1.)*dt*dt*E*y0;
            y = A_c.solve(b);

            b = O*v0 - dt*( E*y0 + E*y )/2.;
            v = O_c.solve(b);

            if(it % 100 == 0)
            {
                std::cout << it << " , jestem w " << t << "\n";
                file1 << t << "\t" << y.transpose()*O*c2 << "\n";
                file2 << t << "\t" << y.transpose()*O*c3 << "\n";
                file3 << t << "\t" << y.transpose()*O*y << "\n";
                file4 << t << "\t" << y.transpose()*E*y << "\n";

                int k = it/100;
                file_a.open(name+std::to_string(k)+".txt");
                computeU(file_a, y);
                file_a.close();
            }

            t += dt; it++;
            y0 = y;
            v0 = v;
        }
        file1.close(); file2.close(); file3.close(); file4.close();
    }

    void computeU(std::ofstream &file, int k);
    void computeU(std::ofstream &file, Eigen::VectorXd y);
    void printMatrix(Eigen::MatrixXd A, int n_rows, int n_cols);
    void printVector(Eigen::VectorXd v, int n_rows);
    void writeMesh(std::ofstream &file);
private:
    int N;

    Node min, max;
    double dx, dy;
    std::vector<Node> nodes;
    std::set<Node> boundaryNodes;
    std::vector<Element> elements;

    std::vector<double> u;

    Eigen::MatrixXd E, O;
    Eigen::VectorXd c;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;

    Eigen::LLT<Eigen::MatrixXd> A_c, O_c;

    void createMesh(int nx);
    std::vector<int> findInnerNodes();

    void initialize();
    void setBoundaries();

};


#endif //LAB04_EQUATION_H
