#include "Grad.hpp"

Point Grad::operator() (const Point& x){  //overwrite of () to switch between the 2 modalities
    if constexpr (Grad_mode == 0){
        Point res(DFex[0](x),DFex[1](x));  //if it's given the exact gradient returns its evaluetions
        return res;
    }
    else{
        return computeGradNum(x);  //else call the function that computes it with finite difference method
    }
}

Point Grad::computeGradNum(const Point& x){  //computes the gradient numerically
    Point hx(h,0.0), hy(0.0,h);  //initialization of the 2 steps in the 2 directions
    double dfdx=(f(x+hx)-f(x-hx))/(2*h), dfdy=(f(x+hy)-f(x-hy))/(2*h);  
    Point res(dfdx,dfdy);

    return res;
}