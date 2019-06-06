//
// Created by Kristina on 06-Nov-18.
//

#ifndef LAB05_CIRCLE_H
#define LAB05_CIRCLE_H

#include <vector>

#include "Node.h"
#include "../Eigen/Eigenvalues"
#include "../Eigen/Dense"

class CircumCircle {
public:
    CircumCircle(Node n1, Node n2, Node n3)
    {
        std::vector<Node> nod{n1, n2, n3};
        Eigen::Matrix3d a, bx, by, c;
        for (int i = 0; i < 3; ++i) {
            bx(i,0) = by(i,0) = c(i,0) = nod[i].x*nod[i].x + nod[i].y*nod[i].y;
            a(i,0) =           by(i,1) = c(i,1) = nod[i].x;
            a(i,1) = bx(i,1)           = c(i,2) = nod[i].y;
            a(i,2) = bx(i,2) = by(i,2)          = 1.;
        }
        double a_det = a.determinant(), bx_det = - bx.determinant(),
                by_det = by.determinant(), c_det = - c.determinant();

        x0 = -bx_det/a_det/2.;
        y0 = -by_det/a_det/2.;
        R = sqrt( bx_det*bx_det + by_det*by_det - 4*a_det*c_det )/fabs(a_det)/2.;
    }

    bool isInside(Node n)
    {
        return pow(n.x - x0, 2) + pow(n.y - y0, 2) < R * R;
    }


    double x0;
    double y0;
    double R;

};


#endif //LAB05_CIRCLE_H
