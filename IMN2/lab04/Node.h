//
// Created by Kristina on 30-Oct-18.
//

#ifndef LAB04_NODE_H
#define LAB04_NODE_H


class Node {
public:
    explicit Node(double xx=0., double yy=0., int n=0)
    : x(xx), y(yy), global_n(n){}

    bool operator<(const Node& n) const
    {
        return global_n < n.global_n;
    }

    double x;
    double y;
    int global_n;
};


#endif //LAB04_NODE_H
