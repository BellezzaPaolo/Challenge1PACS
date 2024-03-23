#include "Parameters.hpp"
#include <fstream>
#include <string>

int main(){
    Point x0(0.0,0.0);  //inital guess of the minimum
    std::vector<std::function<double(Point)>> DFex;    //exact vecctor gradient
    DFex.push_back([](Point x){return x.getY()+16*std::pow(x.getX(),3)+3;});
    DFex.push_back([](Point x){return x.getX()+2*x.getY();});
    double h=1e-6;    //step length for finite differnce method
    std::function<double(Point)> f=[](Point x){return x.getX()*x.getY()+4*std::pow(x.getX(),4)+std::pow(x.getY(),2)+3*x.getX();}; //the main function
    Grad G(DFex,h,f);  //initailization of my grad class
    const Parameters p(x0,1e-6,1e-6,1e6,0.1,0.5,0.2,f,G); //initialization of my parameters class 
    
    p.print();  //show the parameters

    Point res=Minimum(p); //compute the minimum

    std::cout<<"\n\nPoint of minimum: ";
    res.print();
    std::cout<<"With value: "<<f(res)<<std::endl;  //show the results

    return 0;
}