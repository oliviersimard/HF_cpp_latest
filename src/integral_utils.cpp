#include "integral_utils.h"

const int Integrals::_maxLevel = 50;
const int Integrals::_minLevel = 5;

double Integrals::coarse_app(std::function< double(double) > f, double a, double b){
    return (b-a) * (f(a)-f(b))/2.0;
}

double Integrals::trap_app(std::function< double(double) > f, double a, double b){
    double m = (a+b)/2.0;
    return (b-a)/4.0 * (f(a)+2.0*f(m)+f(b));
}

double Integrals::simps_app(std::function< double(double) > f, double a, double b){
    double dx = (b-a)/2.0;
    double m = (b+a)/2.0;
    return dx/3.0 * (f(a)+4.0*f(m)+f(b));
}

double Integrals::trap(std::function< double(double) > f, double a, double b, double tol, int currentlevel){
    double q = trap_app(f,a,b);
    double r = coarse_app(f,a,b);
    if ( (currentlevel>=_minLevel) && (std::abs(q-r)<=1.0*tol) ){
        return q;
    }
    else if (currentlevel>=_maxLevel) {
        std::cout << "Busted the maximum number of iterations allowed!" << std::endl;
        return q;
    }
    else{
        ++currentlevel;
        return ( trap(f,a,(a+b)/2.0,tol,currentlevel) + trap(f,(a+b)/2.0,b,tol,currentlevel) );
    }
}

double Integrals::simps(std::function< double(double) > f, double a, double b, double tol, int currentlevel){
    double q = trap_app(f,a,b);
    double r = coarse_app(f,a,b);
    if ( (currentlevel >= _minLevel) && (std::abs(q-r)<=1.0*tol) ){
        return q;
    }
    else if (currentlevel >= _maxLevel) {
        std::cout << "Busted the maximum number of iterations allowed!" << q << std::endl;
        return q;
    }
    else{
        ++currentlevel;
        return ( simps(f,a,(a+b)/2.0,tol,currentlevel) + simps(f,(a+b)/2.0,b,tol,currentlevel) );
    }
}

double Integrals::integrate_trap(std::function< double(double) > f, double a, double b, double tol){
    return trap(f,a,b,tol,1); // currentlevel starts with number 1, without much surprise.
}

double Integrals::integrate_simps(std::function< double(double) > f, double a, double b, double tol){
    return simps(f,a,b,tol,1);
}

