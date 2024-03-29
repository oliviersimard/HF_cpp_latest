#include "green_utils.h"

using namespace Hubbard;

const std::complex<double> im(0.0,1.0);
const arma::Mat< std::complex<double> > II_(2, 2, arma::fill::eye);
const arma::Mat< std::complex<double> > ZEROS_(2, 2, arma::fill::zeros);

std::ostream& operator<<(std::ostream& os, const FunctorBuildGk& obj){
    return os << "The content is: \n" << "U: " << obj._u << std::endl << "mu: " << obj._mu << std::endl <<
    "beta: " << obj._beta << std::endl << "n_do: " << obj._ndo << std::endl <<
    "size Matsubara arr: " << obj._size << std::endl << "gridK: " << obj._Nk << std::endl;
}

FunctorBuildGk::FunctorBuildGk(double mu,int beta,double u,double ndo,std::vector<double> kArr,std::vector<double> kArr_l,int Nit,int Nk,std::vector< std::complex<double> >& Gup_k) : 
_mu(mu), _u(u), _ndo(ndo), _beta(beta), _Nit(Nit), _Nk(Nk), _kArr(kArr), _kArr_l(kArr_l){
    this->_Gup_k = &Gup_k.front(); 
    this->_size = Gup_k.size();
    for (int j=0; j<_size; j++){
        this->_precomp_wn.push_back(w(j,0.0,_beta));
        this->_precomp_qn.push_back(q(j,_beta));
    }
}

FunctorBuildGk::FunctorBuildGk(MembCarrier* MembCarrierObj,double mu,double ndo,std::vector<double> kArr,std::vector<double> kArr_l,std::vector< std::complex<double> >& Gup_k) : 
_mu(mu), _ndo(ndo), _kArr(kArr), _kArr_l(kArr_l){
    // Initialization from MembCarrier object.
    this->_u=(MembCarrierObj->db_ptr)[0];
    this->_beta=(MembCarrierObj->db_ptr)[2];
    this->_Nit=(MembCarrierObj->int_ptr)[1];
    this->_Nk=(MembCarrierObj->int_ptr)[3];
    // Basically same thing as other constructor.
    this->_Gup_k = &Gup_k.front(); 
    this->_size = Gup_k.size();
    for (int j=0; j<_size; j++){
        this->_precomp_wn.push_back(w(j,0.0,_beta));
        this->_precomp_qn.push_back(q(j,_beta));
    }
    
}

arma::Mat< std::complex<double> > FunctorBuildGk::operator()(int j, double kx, double ky){
    return buildGkAA_2D(j,_mu,_beta,_u,_ndo,kx,ky);
}

arma::Mat< std::complex<double> > FunctorBuildGk::operator()(std::complex<double> w, double kx, double ky){
    return buildGkAA_2D_w(w,_mu,_u,_ndo,kx,ky);
}

inline double FunctorBuildGk::epsk2D(double kx, double ky){
    return -2.0*(cos(kx)+cos(ky));
}

inline double FunctorBuildGk::epsk1D(double kx){
    return -2.0*cos(kx);
}

inline std::complex<double> FunctorBuildGk::w(int n, double mu, int beta){
    return std::complex<double>(0.0,1.0)*(2.0*(double)n+1.0)*M_PI/(double)beta + mu;
}

inline std::complex<double> FunctorBuildGk::q(int n, int beta){
    return std::complex<double>(0.0,1.0)*(2.0*(double)n)*M_PI/(double)beta;
}

arma::Mat< std::complex<double> > FunctorBuildGk::buildGkAA_2D(int j, double mu, int beta, double u, double ndo, double kx, double ky){
    statMat(0,0) = 1.0/( w(j,mu,beta) - u*ndo - epsk2D(kx,ky)*epsk2D(kx,ky)/( w(j,mu,beta) - u*(1.0-ndo) ) ); // G^{AA}_{up} or G^{BB}_{down}
    statMat(1,1) = 1.0/( w(j,mu,beta) - u*(1.0-ndo) - epsk2D(kx,ky)*epsk2D(kx,ky)/( w(j,mu,beta) - u*(ndo) ) ); // G^{AA}_{down} or G^{BB}_{up}
    statMat(0,1) = 0.0+0.0*im; statMat(1,0) = 0.0+0.0*im;
    return statMat;
}

arma::Mat< std::complex<double> > FunctorBuildGk::buildGkAA_2D_w(std::complex<double> w, double mu, double u, double ndo, double kx, double ky){
    // statMat(0,0) = 1.0/( w + mu - u*ndo - epsk2D(kx,ky)*epsk2D(kx,ky)/( w + mu - u*(1.0-ndo) ) ); // G^{AA}_{up} or G^{BB}_{down}
    // statMat(1,1) = 1.0/( w + mu - u*(1.0-ndo) - epsk2D(kx,ky)*epsk2D(kx,ky)/( w + mu - u*(ndo) ) ); // G^{AA}_{down} or G^{BB}_{up}
    statMat(0,0) = 1.0/( w + mu - u*ndo - epsk2D(kx,ky) ); // G^{AA}_{up} or G^{BB}_{down}
    statMat(1,1) = 1.0/( w + mu - u*(1.0-ndo) - epsk2D(kx,ky) ); // G^{AA}_{down} or G^{BB}_{up}
    statMat(0,1) = 0.0+0.0*im; statMat(1,0) = 0.0+0.0*im;
    return statMat;
}

arma::Mat< std::complex<double> > FunctorBuildGk::buildGkAA_1D(int j, double mu, int beta, double u, double ndo, double kx){
    statMat(0,0) = 1.0/( w(j,mu,beta) - u*ndo - epsk1D(kx)*epsk1D(kx)/( w(j,mu,beta) - u*(1.0-ndo) ) ); // G^{AA}_{up}
    statMat(1,1) = 1.0/( w(j,mu,beta) - u*(1.0-ndo) - epsk1D(kx)*epsk1D(kx)/( w(j,mu,beta) - u*(ndo) ) ); // G^{AA}_{down}
    statMat(0,1) = 0.0+0.0*im; statMat(1,0) = 0.0+0.0*im;
    return statMat;
}

arma::Mat< std::complex<double> > FunctorBuildGk::buildGkAA_1D_w(std::complex<double> w, double mu, double u, double ndo, double kx){
    // statMat(0,0) = 1.0/( w + mu - u*ndo - epsk1D(kx)*epsk1D(kx)/( w + mu - u*(1.0-ndo) ) ); // G^{AA}_{up}
    // statMat(1,1) = 1.0/( w + mu - u*(1.0-ndo) - epsk1D(kx)*epsk1D(kx)/( w + mu - u*(ndo) ) ); // G^{AA}_{down}
    statMat(0,0) = 1.0/( w + mu - epsk1D(kx) - u*ndo );
    statMat(1,1) = 1.0/( w + mu - epsk1D(kx) - u*(1.0-ndo) );
    statMat(0,1) = 0.0+0.0*im; statMat(1,0) = 0.0+0.0*im;
    return statMat;
}

arma::Mat< std::complex<double> > FunctorBuildGk::buildGkBB_1D(int j, double mu, int beta, double u, double ndo, double kx){
    arma::Mat< std::complex<double> > tmp_mat(2,2);
    tmp_mat = buildGkAA_1D(j,mu,beta,u,ndo,kx);
    swap(tmp_mat);
    return tmp_mat;
}

arma::Mat< std::complex<double> > FunctorBuildGk::buildGkBB_1D_w(std::complex<double> w, double mu, double u,double ndo, double kx){
    arma::Mat< std::complex<double> > tmp_mat(2,2);
    tmp_mat = buildGkAA_1D_w(w, mu, u, ndo, kx);
    swap(tmp_mat);
    return tmp_mat;
}

arma::Mat< std::complex<double> >& FunctorBuildGk::swap(arma::Mat< std::complex<double> >& M){
    arma::Mat< std::complex<double> > tmp_mat = ZEROS_;
    tmp_mat(0,0) = M(1,1);
    tmp_mat(1,1) = M(0,0);
    M = tmp_mat;
    return M;
}

void FunctorBuildGk::get_ndo_1D(){
    for (int i=0; i<_Nit; i++) {
            
        double ndo_av=0.0;
        for (int kkx=0; kkx<=_Nk; kkx++) {
    
            // calculate Gup_k in Matsubara space (AFM)
            for (int jj=0; jj<_size; jj++){
                *(_Gup_k+jj) = buildGkAA_1D(jj,_mu,_beta,_u,_ndo,_kArr_l[kkx])(0,0);
                //*(_Gup_k+jj) = buildGkBB_1D(jj,_mu,_beta,_u,_ndo,_kArr[kkx])(0,0); // Modified here to BB
            }

            // calculate ndo_k
            double ndo_k=0;
            for (int jj=0; jj<_size; jj++)
                ndo_k += (2./_beta)*( *(_Gup_k+jj)-1./w(jj,0.0,_beta) ).real();
            ndo_k -= 0.5;
            ndo_k *= (-1);
            
            if ((kkx==0) || (kkx==_Nk)){
                ndo_av += 0.5*ndo_k;
            }
            else{
                ndo_av += ndo_k;
            }
            
        }
        ndo_av /= (_Nk);
        _ndo = ndo_av;
    }
}

void FunctorBuildGk::get_ndo_2D(){
    for (int i=0; i<_Nit; i++) {
        
        double ndo_av=0.0;
        for (int kkx=0; kkx<=_Nk; kkx++) {

            for (int kky=0; kky<=_Nk; kky++){
    
                // calculate Gup_k in Matsubara space (AFM)
                for (int jj=0; jj<_size; jj++)
                    *(_Gup_k+jj) = buildGkAA_2D(jj,_mu,_beta,_u,_ndo,_kArr[kkx],_kArr[kky])(0,0);
                // calculate ndo_k
                double ndo_k=0;
                for (int jj=0; jj<_size; jj++)
                    ndo_k += (2./_beta)*( *(_Gup_k+jj)-1./w(jj,0.0,_beta) ).real();
                ndo_k -= 0.5;
                ndo_k *= (-1);
            
                if ( (kky==0) || (kky==_Nk) || (kkx==0) || (kkx==_Nk) ){
                    if ( ((kkx==0) || (kkx==_Nk)) && ((kky==0) || (kky==_Nk)) ){
                        ndo_av += 0.25*ndo_k;
                    }
                    ndo_av += 0.5*ndo_k;
                }
                else
                    ndo_av += ndo_k;
            
            }
        }
        ndo_av /= (_Nk*_Nk);
        _ndo = ndo_av;
    }
}

double FunctorBuildGk::get_double_occupancy_AA(){
/* This function computes U<n_{up}n_{down}> = 1/(\beta*V)*\sum_{k,ikn} \Sigma(k,ikn)G(k,ikn)e^{-ikn0^-} */
    double d=0.0;
    for (int j=0; j<_size; j++){
        std::complex<double> d_k(0.0,0.0);
        std::complex<double> wj = w(j,_mu,_beta);
        for (int k=0; k<_kArr.size(); k++){ // Summing over G(k)
            double epsk = -2.0*cos(_kArr[k]);
            if ( (k==0) || (k==_Nk) ){
                d_k += 0.5*_u*_ndo/( wj - epsk - _u*_ndo );
            }
            else{
                d_k += 1.0*_u*_ndo/( wj - epsk - _u*_ndo);
            }
        }
        d_k /= (_Nk);
        d += (2.0/_beta)*( d_k ).real();
    }
    d *= 1.0/_u;

    return d;
}

K_1D K_1D::operator+(const K_1D& rhs) const{
    K_1D obj(_qx,_iwn);
    obj._qx = this->_qx + rhs._qx;
    obj._iwn = this->_iwn + rhs._iwn;

    return obj;
}

K_1D K_1D::operator-(const K_1D& rhs) const{
    K_1D obj(_qx,_iwn);
    obj._qx = this->_qx - rhs._qx;
    obj._iwn = this->_iwn - rhs._iwn;

    return obj;
}

K_2D K_2D::operator+(const K_2D& rhs) const{
    K_2D obj(_qx,_qy,_iwn);
    obj._qx = this->_qx + rhs._qx;
    obj._qy = this->_qy + rhs._qy;
    obj._iwn = this->_iwn + rhs._iwn;

    return obj;
}

K_2D K_2D::operator-(const K_2D& rhs) const{
    K_2D obj(_qx,_qy,_iwn);
    obj._qx = this->_qx - rhs._qx;
    obj._qy = this->_qy - rhs._qy;
    obj._iwn = this->_iwn - rhs._iwn;

    return obj;
}
