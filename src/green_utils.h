#ifndef GREEN_UTILS_H_
#define GREEN_UTILS_H_

#include<iostream>
#include<complex>
#include<vector>
#include<armadillo>
#include "json_utils.h"

#define SPINDEG 2.0

extern const std::complex<double> im; // Maybe avoid so many const declaration to speed up code.
extern const arma::Mat< std::complex<double> > II_;
extern const arma::Mat< std::complex<double> > ZEROS_;
static arma::Mat< std::complex<double> > statMat(2,2);

class Susceptibility; // class Susceptibility is member of the global scope.
class FFT; // class FFT is member of the global scope.
namespace ThreadFunctor { class ThreadWrapper; }
namespace ThreadFunctor { class ThreadFunctor1D; }
namespace Hubbard { class FunctorBuildGk; } // Have to forward declare class in namespace to be able to overload operator<<.
std::ostream& operator<<(std::ostream&, const Hubbard::FunctorBuildGk&);
void saveGF_grid(const std::string, Hubbard::FunctorBuildGk&);

namespace Hubbard{

class FunctorBuildGk{
    friend std::ostream& ::operator<<(std::ostream& os, const FunctorBuildGk& obj);
    friend void ::saveGF_grid(std::string filename, FunctorBuildGk& obj);
    friend class ::Susceptibility;
    friend class ::FFT;
    friend class ::ThreadFunctor::ThreadFunctor1D;
    friend class ::ThreadFunctor::ThreadWrapper;
    public:
        FunctorBuildGk(double,int,double,double,std::vector<double>,std::vector<double>,int,int,std::vector< std::complex<double> >&);
        FunctorBuildGk(::MembCarrier* MemObj,double mu,double ndo,std::vector<double> kArr,std::vector<double> kArr_l,std::vector< std::complex<double> >& Gup_k);
        FunctorBuildGk()=default;
        ~FunctorBuildGk()=default;
        
        arma::Mat< std::complex<double> > operator()(int, double, double);
        void get_ndo_2D();
        arma::Mat< std::complex<double> > operator()(std::complex<double>, double, double);
        inline arma::Mat< std::complex<double> > operator()(int j, double kx){
            return buildGkAA_1D(j,_mu,_beta,_u,_ndo,kx);
        } // This inline function is called in other translation units!
        void get_ndo_1D();
        inline arma::Mat< std::complex<double> > operator()(std::complex<double> w, double kx){
            return buildGkAA_1D_w(w,_mu,_u,_ndo,kx);
        } // This inline functions are called in other translation units!

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

        double get_double_occupancy_AA();
        inline double get_ndo(){
            return this->_ndo;
        } // Called in main.

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
    K_1D()=default;
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
    K_2D()=default;
    ~K_2D()=default;
    K_2D operator+(const K_2D& rhs) const;
    K_2D operator-(const K_2D& rhs) const;
    
    double _qx, _qy;
    std::complex<double> _iwn;
};

} /* Hubbard */

#endif /* GREEN_UTILS_H_ */