#include "Point.hpp"

Point operator* (double h,const Point& a){  // computes the moltiplication element by element of a double and a point
    Point res(h*a.getX(),h*a.getY());
    return res;
}


Point operator+ (const Point& a,const Point b){   //computes the sum of 2 points without changing their value and returns the result
    Point res(a.getX()+b.getX(),a.getY()+b.getY());
    return res;
}

Point operator- (const Point& a,const Point& b){   //computes the difference of 2 points without changing their value and returns the result
    Point res(a.getX()-b.getX(),a.getY()-b.getY());
    return res;
}

double norm(const Point& a){  //computes the euclidean norm of a point
    return std::sqrt(std::pow(a.getX(),2)+std::pow(a.getY(),2));
}
