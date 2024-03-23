#include <iostream>
#include <cmath>
#ifndef Point_library
#define Point_library

class Point{         //definition of the class that implements a point in a 2d space
    private:
        double x;
        double y;
    
    public:
        Point(double X,double Y): x(X),y(Y){};     //its constructor
        
        void setX(double X){x=X; return;};         //its setter to set to the private members
        void setY(double Y){y=Y; return;};

        double getX()const{return x;};          //its getter to acces the private members
        double getY()const{return y;};

        void print() const{                       //a function to print the point
            std::cout<<"( "<<x<<" , "<<y<<" )"<<std::endl;
            return;
        }

};

//overwrite of all the needed operator
Point operator* (double h,const Point& a);
Point operator+ (const Point& a,const Point b);
Point operator- (const Point& a,const Point& b);
double norm(const Point& a);

#endif