//
// Created by Kristina on 30-Oct-18.
//

#include "Element.h"

Element::Element(int i, int j, int nx, double dx, double dy, Node _min, Node _max)
                : min_global(_min), max_global(_max)
{

    nodes.emplace_back(dx*(i-1), dy*(j-1), i+(j-1)*nx); // l=1
    nodes.emplace_back(dx*i, dy*(j-1), (i+1)+(j-1)*nx); // l=2
    nodes.emplace_back(dx*i, dy*j, (i+1)+j*nx); // l=3
    nodes.emplace_back(dx*(i-1), dy*j, i+j*nx); // l=4

    min = nodes[0], max = nodes[2];
}

Node &Element::get(int l) {
    return nodes[l-1];
}

int Element::get_alpha(int l) {
    if(l==1) return 1;
    if(l==2) return 2;
    if(l==3) return 2;
    if(l==4) return 1;
}

int Element::get_beta(int l) {
    if(l==1) return 1;
    if(l==2) return 1;
    if(l==3) return 2;
    if(l==4) return 2;
}

double Element::fi(int a, int b, double ksi) {
    if(a==0 && b==1) return 1./2. - 3.*ksi/4. + pow(ksi, 3)/4.;
    if(a==0 && b==2) return 1./2. + 3.*ksi/4. - pow(ksi, 3)/4.;
    if(a==1 && b==1) return (1 - ksi - ksi*ksi + pow(ksi, 3)) / 4.;
    if(a==1 && b==2) return (- 1 - ksi + ksi*ksi + pow(ksi, 3)) / 4.;
}

double Element::d_fi(int a, int b, double ksi) {
    double d_ksi = 0.001;
    return ( fi(a,b,ksi+d_ksi) + fi(a,b,ksi-d_ksi) - 2.*fi(a,b,ksi) ) /d_ksi/d_ksi;
}

double Element::Jacobian() {
    return (get(2).x - get(1).x) * (get(4).y - get(1).y) / 4.;
}

double Element::w(int i, double ksi1, double ksi2) {
    if(i==1) return (1 - ksi1)*(1 - ksi2)/4.;
    if(i==2) return (1 + ksi1)*(1 - ksi2)/4.;
    if(i==3) return (1 + ksi1)*(1 + ksi2)/4.;
    if(i==4) return (1 - ksi1)*(1 + ksi2)/4.;
}

double Element::x(double ksi1, double ksi2) {
    double result = 0.;
    for (int i = 1; i <= 4; ++i) {
        result += get(i).x * w(i,ksi1,ksi2);
    }
    return result;
}

double Element::y(double ksi1, double ksi2) {
    double result = 0.;
    for (int i = 1; i <= 4; ++i) {
        result += get(i).y * w(i,ksi1,ksi2);
    }
    return result;
}

double Element::integrate(int l1, int i1, int i2, double (*ro)(double,double)) {
    std::vector<double> weights = gl5.get_weights();
    std::vector<double> ksi = gl5.get_roots();
    double result = 0.;

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++) {
            result += +weights[i] * weights[j] * fi(i1, get_alpha(l1), ksi[i]) * fi(i2, get_beta(l1), ksi[j])
                        * Jacobian() * ro( x(ksi[i],ksi[j]) , y(ksi[i],ksi[j]) ) * (-1.);
        }
    }
    return result;
}

double Element::integrate(int l1, int l2, int i1, int i2, int j1, int j2) {
    std::vector<double> weights = gl5.get_weights();
    std::vector<double> ksi = gl5.get_roots();
    double result = 0.;

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            result += weights[i] * weights[j] * fi(i1,get_alpha(l1),ksi[i]) * fi(i2,get_beta(l1),ksi[j]) *
                      ( fi(j2,get_beta(l2),ksi[j]) * d_fi(j1,get_alpha(l2),ksi[i]) +
                        fi(j1,get_alpha(l2),ksi[i]) * d_fi(j2,get_beta(l2),ksi[j]) );
        }
    }
    return result;
}

void Element::matrixForNodes(int l1, int l2, Eigen::MatrixXd &S_node) {
    int i, j;
    for (int k = 0; k < 4; ++k) {
        for (int l = 0; l < 4; ++l) {
            S_node(k,l) = 0.;
        }
    }

    for(int i1=0; i1<=1; i1++) {
        for(int j1=0; j1<=1; j1++) {
            for(int i2=0; i2<=1; i2++) {
                for(int j2=0; j2<=1; j2++)
                {
                    i = i1 + 2*i2;
                    j = j1 + 2*j2;
                    S_node(i, j) = integrate(l1, l2, i1, i2, j1, j2);
                }
            }
        }
    }
}

void Element::addToGlobalMatrix(Eigen::MatrixXd &S) {
    int p0, q0;
    Eigen::MatrixXd S_node(4,4);
    for(int l1=1; l1<=4; l1++) {
        for(int l2=1; l2<=4; l2++)
        {
            p0 = 4*(get(l1).global_n - 1);
            q0 = 4*(get(l2).global_n - 1);
            //std::cout << "gloab_n = " << get(l1).global_n << "p0 = " << p0 << "\n";
            matrixForNodes(l1, l2, S_node);
            // merging global S and S for_nodes
            for(int i=0; i<4; i++) {
                for(int j=0; j<4; j++) {
                    S(p0+i, q0+j) += S_node(i,j);
                }
            }
        }
    }
}

void Element::addToGlobalVector(Eigen::VectorXd &F, double (*ro)(double, double)) {
    int p;
    for(int l1=1; l1<=4; l1++) {
        for(int i1=0; i1<=1; i1++)
        {
            for(int i2=0; i2<=1; i2++) {
                p = 4 * (get(l1).global_n - 1) + i1 + 2*i2;
                //std::cout << "p = " << p << "\n";
                F(p) += integrate(l1, i1, i2, ro);
            }
        }
    }
}

std::vector<Node> Element::getEdgeNodes() {
    std::vector<Node> edgeNodes;
    for (int i = 1; i <= 4; ++i) {
        if( get(i).x == max_global.x || get(i).x == min_global.x || get(i).y == max_global.y || get(i).y == min_global.y)
        {
            edgeNodes.push_back(get(i));
        }
    }
    return std::move(edgeNodes);
}

bool Element::isElement(double x, double y) {
    return (x >= min.x && x < max.x) && (y >= min.y && y < max.y);
}

double Element::ksi1(double x) {
    return 2. * ( x - (get(2).x + get(1).x)/2. ) / (get(2).x - get(1).x);
}

double Element::ksi2(double y) {
    return 2. * ( y - (get(4).y + get(1).y)/2. ) / (get(4).y - get(1).y);
}

double Element::computeUValue(double x, double y, Eigen::VectorXd &c) {
    double result = 0.;
    int p;
    for(int l=1; l<=4; l++)
    {
        for (int i1 = 0; i1 <= 1; ++i1) {
            for (int i2 = 0; i2 <= 1 ; ++i2) {
                p = 4*( get(l).global_n - 1 ) + i1 + 2*i2;
                result += c(p) * fi(i1, get_alpha(l), ksi1(x)) * fi(i2, get_beta(l), ksi2(y));
            }
        }
    }
    return result;
}





