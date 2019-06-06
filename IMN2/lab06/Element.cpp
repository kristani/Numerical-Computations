//
// Created by Kristina on 30-Oct-18.
//

#include "Element.h"

Element::Element(Node *n1, Node *n2, Node *n3, Node *_min, Node *_max)
        : min_global(_min), max_global(_max), H_local(3,3), M_local(3,3)
{
    nodes.push_back(n1);
    nodes.push_back(n2);
    nodes.push_back(n3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            H_local(i,j) = 0.;
            M_local(i,j) = 0.;
        }
    }
    //std::cout << "J = "<< jacobian << "\n";
}


Node& Element::get(int l) {
    return *nodes[l-1];
}

double Element::computeArea(Node n1, Node n2, Node n3) {
    //double a = -n1.x*(n3.y-n2.y)+x2*y3-x3*y2+(x3-x2)*y1;
    double a = (n2.x - n1.x)*(n3.y - n1.y) - (n3.x - n1.x)*(n2.y - n1.y);
    return a/2.;
}

bool Element::isInside(double xx, double yy) {
    Node n(xx,yy);
    double area1 = fabs( computeArea(n, *nodes[1], *nodes[2]) );
    double area2 = fabs( computeArea(*nodes[0], n, *nodes[2]) );
    double area3 = fabs( computeArea(*nodes[0], *nodes[1], n) );
    double a = area1 + area2 + area3;
    return fabs((area1 + area2 + area3) - area) < 1e-8 ;
}

double Element::fi(int i, double ksi1, double ksi2) {
    if(i==0) return -(ksi2 + ksi1)/2.;
    if(i==1) return (1 + ksi1)/2.;
    if(i==2) return (ksi2 + 1)/2.;
}

double Element::d_fi(int flag, int i, double ksi1, double ksi2) // flag=1 - ksi1, flag=2 - ksi2
{
    double d_ksi = 1e-4;
    if(flag == 1) // d_fi po d_ksi1
    { return ( fi(i,ksi1+d_ksi,ksi2) - fi(i,ksi1-d_ksi,ksi2) ) /d_ksi/2.; }

    if(flag == 2) // d_fi po d_ksi2
    { return ( fi(i,ksi1,ksi2+d_ksi) - fi(i,ksi1,ksi2-d_ksi) ) /d_ksi/2.; }
}

double Element::x(double ksi1, double ksi2) {
    double result = 0.;
    for (int i = 0; i <= 2; ++i) {
        result += nodes[i]->x * fi(i,ksi1,ksi2);
    }
    return result;
}

double Element::y(double ksi1, double ksi2) {
    double result = 0.;
    for (int i = 0; i <= 2; ++i) {
        result += nodes[i]->y * fi(i,ksi1,ksi2);
    }
    return result;
}

double Element::d_x(int flag, double ksi1, double ksi2) {
    double d_ksi = 1e-4;
    if(flag == 1) // d_x po d_ksi1
    { return ( x(ksi1+d_ksi,ksi2) - x(ksi1-d_ksi,ksi2) ) /d_ksi/2.; }

    if(flag == 2) // d_x po d_ksi2
    { return ( x(ksi1,ksi2+d_ksi) - x(ksi1,ksi2-d_ksi) ) /d_ksi/2.; }
}

double Element::d_y(int flag, double ksi1, double ksi2) {
    double d_ksi = 1e-4;
    if(flag == 1) // d_x po d_ksi1
    { return ( y(ksi1+d_ksi,ksi2) - y(ksi1-d_ksi,ksi2) ) /d_ksi/2.; }

    if(flag == 2) // d_x po d_ksi2
    { return ( y(ksi1,ksi2+d_ksi) - y(ksi1,ksi2-d_ksi) ) /d_ksi/2.; }
}

double Element::d_ksi1(int flag, double ksi1, double ksi2) {
    // d_ksi1/dx = dy / d_ksi2
    if(flag == 1)   return d_y(2,ksi1, ksi2)/jacobian;
    // d_ksi1/dy = -dx / d_ksi2
    if(flag == 2)   return - d_x(2,ksi1, ksi2)/jacobian;
}

double Element::d_ksi2(int flag, double ksi1, double ksi2) {
    // d_ksi2/dx = -dy/dksi_1
    if(flag == 1)   return - d_y(1,ksi1, ksi2)/jacobian;
    // d_ksi2/dy = dx/dksi_1
    if(flag == 2)   return d_x(1,ksi1, ksi2)/jacobian;
}

double Element::Jacobian() {
    double ksi1 = -0.5, ksi2 = -0.5;
    double dx1 = d_x(1,ksi1,ksi2);
    double dy2 = d_y(2,ksi1,ksi2);
    double dx2 = d_x(2,ksi1,ksi2);
    double dy1 = d_y(1,ksi1,ksi2);
    double d = dx1 * dy2 - dy1 * dx2;
//    double r = ((-nodes[0]->x + nodes[1]->x)*(-nodes[0]->y + nodes[2]->y) - (-nodes[0]->y + nodes[1]->y)*(-nodes[0]->x + nodes[2]->x))/4.;
//    double j = fabs(-nodes[0]->x*(nodes[2]->y-nodes[1]->y) + nodes[1]->x*nodes[2]->y - nodes[2]->x*nodes[1]->y + (nodes[2]->x-nodes[1]->x)*nodes[0]->y )*0.5/2.;
    return d;
}

double Element::nabla_fi(int i, int j, double ksi1, double ksi2) {
    std::vector<double> fi_i(2);
    fi_i[0] = d_fi(1,i,ksi1,ksi2)*d_ksi1(1,ksi1,ksi2) + d_fi(2,i,ksi1,ksi2)*d_ksi2(1,ksi1,ksi2);
    fi_i[1] = d_fi(1,i,ksi1,ksi2)*d_ksi1(2,ksi1,ksi2) + d_fi(2,i,ksi1,ksi2)*d_ksi2(2,ksi1,ksi2);

    std::vector<double> fi_j(2);
    fi_j[0] = d_fi(1,j,ksi1,ksi2)*d_ksi1(1,ksi1,ksi2) + d_fi(2,j,ksi1,ksi2)*d_ksi2(1,ksi1,ksi2);
    fi_j[1] = d_fi(1,j,ksi1,ksi2)*d_ksi1(2,ksi1,ksi2) + d_fi(2,j,ksi1,ksi2)*d_ksi2(2,ksi1,ksi2);

    return fi_i[0]*fi_j[0] + fi_i[1]*fi_j[1];
}

void Element::checkLocalIndexes() {
    double w  = computeArea(*nodes[0],*nodes[1],*nodes[2]);
    if(w < 0) {
        std::swap(nodes[0], nodes[2]);}
}

double Element::integrateH(int i, int j) {
    double ksi1[]={-0.3333333, -0.0597158717,-0.0597158717, -0.8805682564, -0.7974269853, -0.7974269853, 0.5948539707};
    double ksi2[]={-0.3333333, -0.0597158717,-0.8805682564, -0.0597158717, -0.7974269853,  0.5948539707,-0.7974269853};
    double weights[]={0.45,0.2647883055, 0.2647883055, 0.2647883055, 0.251878361, 0.251878361, 0.251878361};

    double result = 0.;

    for(int k = 0; k < 7; k++) {
        result += weights[k] * jacobian * nabla_fi(i,j,ksi1[k],ksi2[k]);
    }
    return result;
}

double Element::integrateM(int i, int j) {
    double ksi1[]={-0.3333333, -0.0597158717,-0.0597158717, -0.8805682564, -0.7974269853, -0.7974269853, 0.5948539707};
    double ksi2[]={-0.3333333, -0.0597158717,-0.8805682564, -0.0597158717, -0.7974269853,  0.5948539707,-0.7974269853};
    double weights[]={0.45,0.2647883055, 0.2647883055, 0.2647883055, 0.251878361, 0.251878361, 0.251878361};

    double result = 0.;

    for(int k = 0; k < 7; k++) {
        result += weights[k] * jacobian * fi(i,ksi1[k],ksi2[k]) * fi(j,ksi1[k],ksi2[k]);
    }
    return result;
}

void Element::fillLocalStiffnessMatrix() {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            H_local(i,j) = integrateH(i, j);
        }
    }

    //std::cout << "\nMacierz E:\n";
    //std::cout << H_local;
}
void Element::fillLocalMassMatrix() {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            M_local(i,j) = integrateM(i, j);
        }
    }

    //std::cout << "\nWektor F:\n";
    //std::cout << F_local;
}

void Element::fillGlobalStiffnessMatrix(Eigen::MatrixXd &H) {
    int p, q;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            p = nodes[i]->global_n - 1;
            q = nodes[j]->global_n - 1;
            H(p,q) += H_local(i,j);
        }
    }
}

void Element::fillGlobalMassMatrix(Eigen::MatrixXd &M) {
    int p, q;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            p = nodes[i]->global_n - 1;
            q = nodes[j]->global_n - 1;
            M(p,q) += M_local(i,j);
        }
    }
}

double Element::ksi2(double x, double y) {
    double x1 = nodes[0]->x, y1 = nodes[0]->y;
    double x2 = nodes[1]->x, y2 = nodes[1]->y;
    double x3 = nodes[2]->x, y3 = nodes[2]->y;
    return -1./(-x1*y3-y1*x2+x2*y3+y1*x3+x1*y2-y2*x3)*(2*x1*y-x1*y2-x1*y3+y1*x3-2*x2*y+x2*y3-y2*x3-2*y1*x+y1*x2+2*y2*x);
}

double Element::ksi1(double x, double y) {
    double x1 = nodes[0]->x, y1 = nodes[0]->y;
    double x2 = nodes[1]->x, y2 = nodes[1]->y;
    double x3 = nodes[2]->x, y3 = nodes[2]->y;
    return (2*x1*y-x1*y2-x1*y3+y1*x2+y2*x3+y1*x3-2*y1*x-2*x3*y+2*x*y3-x2*y3)/(-x1*y3-y1*x2+x2*y3+y1*x3+x1*y2-y2*x3);
}

double Element::computeUValue(double x, double y, Eigen::VectorXd &c) {
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

void Element::print() {
    for (auto &node : nodes) {
        std::cout << "(" << node->x << " , " << node->y << ") " << node->global_n << "\t";
    }
    std::cout << "\n";
}
