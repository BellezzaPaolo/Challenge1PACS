#include "Point.hpp"
#include <vector>
#include <functional>
#define Grad_mode 1
//0 if it's given the exact one, 1 if it's used the finite difference method to compute the gradient
#ifndef Grad_library
#define Grad_library


class Grad{   //definition of the gradient and the methods allow to make it exact or with finite differences
    private:
        std::vector<std::function<double(Point)>> DFex;   //a vector conteining the 2 exact component of the gradient
        double h;                                         //step length
        std::function<double(Point)> f;                   //the function

    public:
        Grad(std::vector<std::function<double(Point)>> DFEX,double H,std::function<double(Point)> F):DFex(DFEX),h(H),f(F){};  //constructor

        Point operator() (const Point& x);  //overwrite of () to switch between the 2 modalities
        Point computeGradNum(const Point& x);  //computes the gradient numerically

};
#endif