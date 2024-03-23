#include "Parameters.hpp"
//#define S 1

bool condition(const Point& x, double alpha,double omega,std::function<double(Point)> f,Grad G){
    return f(x)-f(x-alpha*G(x))>omega*alpha*std::pow(norm(G(x)),2); //check of the condition in Armijo rule and returns true if satisfied
}

template <strategy S>
double calcAlphak(const Parameters& p,double k,const Point& x){   //computes the alphaK at the k-th iteration dipending on the strategy choosen
    double alphaK;
    if constexpr (S==strategy::Inverse)         //if the choosen strategy is inverse decay camputes alphaK accordingly
        alphaK=p.getAlpha0()/(1+p.getMu()*k);
    else if constexpr (S==strategy::Armijo){   //else computes alphaK using Armijo rule
        alphaK=p.getAlpha0();
        while(!condition(x, alphaK,p.getOmega(),p.getF(),p.getG()) && alphaK>1e-3){  //check of the condition of arrest
            alphaK/=2;
        }
    }
    return alphaK;
}


Point Minimum(const Parameters& p){
    Grad G=p.getG();    //gradient of the function
   double alphaK=p.getAlpha0();
    Point x0=p.getX0(),xN=x0-alphaK*G(x0);  //first 2 iterate of the method
    bool condition=true;
    double k=0;  //counter of iterations
 
    const strategy S=strategy::Inverse;  //defines the strategy


    while(condition){
        k++;           //updates the counter

        alphaK=calcAlphak<S>(p,k,x0);    //call the functon to compute alphaK
        xN=x0-alphaK*(G(x0));            //update the guess of the minimum

        if(norm(xN-x0)<p.getTolls() || norm(G(xN))<p.getTollr() || k>p.getMaxIter()){  //check of the stop condition
            condition=false;
        }
        x0=xN;      //update of the previous guess before restart
    }

    std::cout<<"iterate: "<< k<<std::endl;  //shows the number of iterations
    return xN;
}


void Parameters::print() const{   //displaies all the parameters

    std::cout<<"\n##########################################\n##########################################"<<std::endl;
    std::cout<<"the parameters of my problem are:"<<std::endl;
    std::cout<<"-initial point: ";
    x0.print();
    std::cout<<"-tollerance over residual: "<<tollr<<std::endl;
    std::cout<<"-tollerance over gradient: "<<tolls<<std::endl;
    std::cout<<"-maximum number of iterations: "<<maxIter<<std::endl;
    std::cout<<"-initial value of alpha: "<<alpha0<<std::endl;
    std::cout<<"\n##########################################\n##########################################\n"<<std::endl;
}