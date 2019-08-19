#ifndef Integral_utils_H_
#define Integral_utils_H_

#include <armadillo>
#include <complex>
#include <functional>

// #define INTEGRAL

/* Template structure to call functions in classes. */
template<typename T, typename C, typename Q>
struct functorStruct{
    //using matCplx = arma::Mat< std::complex<double> >;
    using funct_init_t = arma::Mat< std::complex<double> > (C::*)(Q model, T kk, T qq, int n, int l);
    using funct_con_t = arma::Mat< std::complex<double> > (C::*)(Q model, T kk, T qq, int n, int l, arma::Mat< std::complex<double> > SE);

    functorStruct(funct_init_t initFunct, funct_con_t conFunct);
    arma::Mat< std::complex<double> > callInitFunct(C& obj, Q model, T kk, T qq, int n, int l);
    arma::Mat< std::complex<double> > callConFunct(C& obj, Q model, T kk, T qq, int n, int l, arma::Mat< std::complex<double> > SE);

    private:
        funct_init_t _initFunct;
        funct_con_t _conFunct;
};

template<typename T, typename C, typename Q>
functorStruct<T,C,Q>::functorStruct(funct_init_t initFunct, funct_con_t conFunct) : _initFunct(initFunct), _conFunct(conFunct){};

template<typename T, typename C, typename Q>
arma::Mat< std::complex<double> > functorStruct<T,C,Q>::callInitFunct(C& obj, Q model, T kk, T qq, int n, int l){
    return (obj.*_initFunct)(model, kk, qq, n, l);
}

template<typename T, typename C, typename Q>
arma::Mat< std::complex<double> > functorStruct<T,C,Q>::callConFunct(C& obj, Q model, T kk, T qq, int n, int l, arma::Mat< std::complex<double> > SE){
    return (obj.*_conFunct)(model, kk, qq, n, l, SE);
}

class Integrals{
    public:
        double coarse_app(std::function< double(double) >,double,double);
        double trap_app(std::function< double(double) >,double,double);
        double simps_app(std::function< double(double) >,double,double);

        double trap(std::function< double(double) >,double,double,double,int);
        double simps(std::function< double(double) >,double,double,double,int);

        double integrate_trap(std::function< double(double) >,double,double,double);
        double integrate_simps(std::function< double(double) >,double,double,double);

    private:
        static const int _maxLevel;
        static const int _minLevel;

};

#endif /* Integral_utils_H_ */