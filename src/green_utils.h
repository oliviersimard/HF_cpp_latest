#ifndef GREEN_UTILS_H_
#define GREEN_UTILS_H_

#include<iostream>
#include<complex>
#include<vector>
#include<armadillo>
#include "json_utils.h"

const std::complex<double> im(0.0,1.0);
const arma::Mat< std::complex<double> > II_(2, 2, arma::fill::eye);
const arma::Mat< std::complex<double> > ZEROS_(2, 2, arma::fill::zeros);
static arma::Mat< std::complex<double> > statMat(2,2);

class Susceptibility; // class Susceptibility is member of the global scope.
class FFT; // class FFT is member of the global scope.
namespace Hubbard { class FunctorBuildGk; } // Have to forward declare class in namespace to be able to overload operator<<.
std::ostream& operator<<(std::ostream& os, const Hubbard::FunctorBuildGk& obj);

namespace Hubbard{

class FunctorBuildGk{
    friend std::ostream& ::operator<<(std::ostream& os, const FunctorBuildGk& obj);
    friend class ::Susceptibility;
    friend class ::FFT;
    public:
        FunctorBuildGk(double,int,double,double,std::vector<double>,std::vector<double>,int,int,std::vector< std::complex<double> >&);
        FunctorBuildGk(::MembCarrier* MemObj,double mu,double ndo,std::vector<double> kArr,std::vector<double> kArr_l,std::vector< std::complex<double> >& Gup_k);
        ~FunctorBuildGk()=default;
        
        arma::Mat< std::complex<double> > operator()(int, double, double);
        void get_ndo_2D();
        arma::Mat< std::complex<double> > operator()(std::complex<double>, double, double);
        arma::Mat< std::complex<double> > operator()(int, double);
        void get_ndo_1D();
        arma::Mat< std::complex<double> > operator()(std::complex<double>, double);

        double epsk2D(double,double);
        double epsk1D(double);
        std::complex<double> w(int,double,int);
        std::complex<double> q(int, int);
        arma::Mat< std::complex<double> >& swap(arma::Mat< std::complex<double> >&);

        arma::Mat< std::complex<double> > buildGkAA_2D(int,double,int,double,double,double,double);
        arma::Mat< std::complex<double> > buildGkAA_2D_w(std::complex<double>,double,double,double,double,double);
        arma::Mat< std::complex<double> > buildGkAA_1D(int,double,int,double,double,double);
        arma::Mat< std::complex<double> > buildGkAA_1D_w(std::complex<double>,double,double,double,double);
        arma::Mat< std::complex<double> > buildGkBB_1D(int,double,int,double,double,double);
        arma::Mat< std::complex<double> > buildGkBB_1D_w(std::complex<double>,double,double,double,double);

        double get_ndo();

    private:
        double _mu, _u, _ndo;
        int _beta, _Nit, _Nk;
        std::vector<double> _kArr, _kArr_l;
        std::complex<double>* _Gup_k;
        size_t _size;
        std::vector< std::complex<double> > _precomp_wn, _precomp_qn;
};

struct K_1D{
    K_1D(double qx, std::complex<double> iwn) : _qx(qx), _iwn(iwn){};
    ~K_1D()=default;
    K_1D operator+(const K_1D& rhs) const;
    K_1D operator-(const K_1D& rhs) const;

    double _qx;
    std::complex<double> _iwn;
};

struct K_2D : K_1D{
    K_2D(double qx, double qy, std::complex<double> iwn) : K_1D(qx,iwn){
        this->_qy = qy;
    }
    ~K_2D()=default;
    K_2D operator+(const K_2D& rhs) const;
    K_2D operator-(const K_2D& rhs) const;
    
    double _qx, _qy;
    std::complex<double> _iwn;
};

} /* Hubbard */

#endif /* GREEN_UTILS_H_ */