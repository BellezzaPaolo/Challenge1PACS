#include "Parameters.hpp"

bool condition(const Point& x, double alpha,double omega,std::function<double(Point)> f,Grad G){
    return f(x)-f(x-alpha*G(x))>omega*alpha*std::pow(norm(G(x)),2); //check of the condition in Armijo rule and returns true if satisfied
}

/*
template <method M>
double calcAlphak(const Parameters& p,double k,const Point& x){   //computes the alphaK at the k-th iteration dipending on the strategy choosen
    double alphaK;
    if constexpr (M==method::heavyBall || M==method::Nesterov)         //if the choosen strategy is inverse decay camputes alphaK accordingly
        alphaK=p.getAlpha0()/(1+p.getMu()*k);
    else if constexpr (M==method::Armijo){   //else computes alphaK using Armijo rule
        alphaK=p.getAlpha0();
        while(!condition(x, alphaK,p.getOmega(),p.getF(),p.getG()) && alphaK>1e-3){  //check of the condition of arrest
            alphaK/=2;
        }
    }
    return alphaK;
}
*/

template <method M>
Point Minimum(const Parameters& p){
    Grad G=p.getG();    //gradient of the function
    double alphaK=p.getAlpha0();
    Point x0=p.getX0(),x1=x0-alphaK*G(x0),x2(1.0,1.0);  //first 2 iterate of the method
    Point y(0.0,0.0); //for nesterov and heavy-ball method
    bool condizione=true;
    double k=0;  //counter of iterations
    //const method M=method::heavyBall;

    while(condizione){
        k++;           //updates the counter

        if constexpr (M==method::Armijo){ //update the guess of the minimum with the right method
            alphaK=p.getAlpha0();
            while(!condition(x1, alphaK,p.getOmega(),p.getF(),p.getG())){  //check of the condition of arrest
                alphaK/=2;
            }

            x2=x1-alphaK*(G(x1));                        
        }
        else if constexpr (M==method::heavyBall){
            alphaK=p.getAlpha0()/(1+p.getMu()*k);

            x2=x1-alphaK*(G(x1))+p.getNu()*(x1-x0);
        }
        else{
            alphaK=p.getAlpha0()/(1+p.getMu()*k);

            y=x1+p.getNu()*(x1-x0);
            x2=y-alphaK*(G(y));
        }
 

        if(norm(x2-x1)<p.getTolls() || norm(G(x2))<p.getTollr() || k>p.getMaxIter()){  //check of the stop condition
            condizione=false;
        }
        x1=x2;
        x0=x1;      //update of the previous guesses before restart
    }

    std::cout<<"iterate: "<< k<<std::endl;  //shows the number of iterations
    return x2;
}

template Point Minimum<method::Armijo>(const Parameters& p);
template Point Minimum<method::heavyBall>(const Parameters& p);
template Point Minimum<method::Nesterov>(const Parameters& p);


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