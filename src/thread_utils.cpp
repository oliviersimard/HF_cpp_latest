#include "thread_utils.h"

using namespace ThreadFunctor;

ThreadFunctor1D::ThreadFunctor1D(std::complex<double>& upper_level,Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D& q,arma::Mat< std::complex<double> >& mat) : 
_upper_level(upper_level), _Gk(Gk), _q(q), _ktb(mat) {}

ThreadFunctor1D::ThreadFunctor1D(ThreadFunctor1D&& rhs) : 
_upper_level(std::move(rhs)._upper_level), _Gk(std::move(rhs)._Gk), _q(std::move(rhs)._q), _ktb(std::move(rhs)._ktb) {
    // Emptying the source object
    rhs.~ThreadFunctor1D();
}

ThreadFunctor1D::ThreadFunctor1D(const ThreadFunctor1D& rhs) : 
 _upper_level(rhs._upper_level), _Gk(rhs._Gk), _q(rhs._q), _ktb(rhs._ktb) {}

void ThreadFunctor1D::operator()(int ktilde, int kbar){
    std::complex<double> tmp_val_kbar=0.0+0.0*im;
    for (int wtilde=0; wtilde<_Gk._size; wtilde++){
        for (int wbar=0; wbar<_Gk._size; wbar++){
            std::complex<double> tmp_val = _Gk(
                    _Gk._precomp_wn[wtilde],_Gk._kArr_l[ktilde]
                    )(0,0)*_Gk(
                    _Gk._precomp_wn[wtilde]-_q._iwn,_Gk._kArr_l[ktilde]-_q._qx
                    )(0,0)*gamma_oneD_spsp( _Gk._kArr_l[ktilde], _Gk._precomp_wn[wtilde], _Gk._kArr_l[kbar], _Gk._precomp_wn[wbar] )*_Gk(
                    _Gk._precomp_wn[wbar]+_q._iwn,_Gk._kArr_l[kbar]+_q._qx
                    )(1,1)*_Gk(
                    _Gk._precomp_wn[wbar],_Gk._kArr_l[kbar]
                    )(1,1);
            tmp_val_kbar += tmp_val;
            std::lock_guard<std::mutex> guard(mutx); // Internally locked
            _upper_level += tmp_val;
        } 
    }
    _ktb(ktilde,kbar)=tmp_val_kbar; 
}

std::complex<double> ThreadFunctor1D::gamma_oneD_spsp(double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar){
    std::complex<double> lower_level=0.0+0.0*im;
    for (int wttilde=0; wttilde<_Gk._size; wttilde++){
        for (size_t qttilde=0; qttilde<_Gk._kArr_l.size(); qttilde++){
            lower_level += _Gk(wtilde-_Gk._precomp_qn[wttilde],ktilde-_Gk._kArr_l[qttilde])(0,0)*_Gk(wbar+_q._iwn-_Gk._precomp_qn[wttilde],kbar+_q._qx-_Gk._kArr_l[qttilde])(1,1);
        }
    }
    lower_level *= -1.0*_Gk._u/(_Gk._beta*_Gk._Nk); /// Removed minus sign
    lower_level += 1.0;
    return _Gk._u/lower_level;
}


void ThreadFunctor1D::join_all(std::vector<std::thread>& grp){
    for (auto& thread : grp){
        if (thread.joinable()){
            thread.join();
        }
    }
}

/* ThreadWrapper member functions */

// ThreadWrapper::ThreadWrapper(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D& q,arma::cx_dmat::iterator matPtr,arma::cx_dmat::iterator matWPtr) : 
// _Gk(Gk), _q(q), _ktb(matPtr), _ktbW(matWPtr) {}

ThreadWrapper::ThreadWrapper(Hubbard::FunctorBuildGk Gk,Hubbard::K_1D& q,double ndo_converged) : _q(q){
    this->_ndo_converged=ndo_converged;
    this->_Gk=Gk;
}

void ThreadWrapper::operator()(int ktilde, int kbar, double beta){
    std::complex<double> tmp_val_kt_kb(0.0,0.0);
    std::complex<double> tmp_val_weigths(0.0,0.0);
    // cout << beta << " " << _Gk._kArr_l[ktilde] << " " << _q._iwn << " " << _q._qx << " " << _Gk._kArr_l[kbar] << endl;
    for (int wtilde=0; wtilde<_Gk._size; wtilde++){
        for (int wbar=0; wbar<_Gk._size; wbar++){
            tmp_val_kt_kb += gamma_oneD_spsp(_Gk._kArr_l[ktilde],std::complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta),_Gk._kArr_l[kbar],std::complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta));
            tmp_val_weigths += buildGK(
                                    std::complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta),_Gk._kArr_l[ktilde]
                                    )[0]*buildGK(
                                    std::complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)-_q._iwn,_Gk._kArr_l[ktilde]-_q._qx
                                    )[0]*buildGK(
                                    std::complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)-_q._iwn,_Gk._kArr_l[kbar]-_q._qx
                                    )[1]*buildGK(
                                    std::complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta),_Gk._kArr_l[kbar]
                                    )[1]; 
            if ((wtilde==0) && (wbar==0)){
                std::cout << "Thread id: " << std::this_thread::get_id() << std::endl;
                std::cout << "ktilde: " << _Gk._kArr_l[ktilde] << std::endl;
                std::cout << "kbar: " << _Gk._kArr_l[kbar] << std::endl;
                std::cout << "gamma_oneD_spsp: " << tmp_val_kt_kb << std::endl;
                std::cout << "weights: " << tmp_val_weigths << std::endl;
            }
        } 
    }
    // lock_guard<mutex> guard(mutx);
    matGamma(kbar,ktilde) = tmp_val_kt_kb; // These matrices are static variables.
    matWeigths(kbar,ktilde) = 1.0/(_Gk._beta*_Gk._Nk)/(_Gk._beta*_Gk._Nk)*tmp_val_weigths;
    //cout << "Gamma for " << "ktilde " << ktilde << " and kbar " << kbar << ": " << matGamma(kbar,ktilde) << "\n";
    //cout << "Weigths for " << "ktilde " << ktilde << " and kbar " << kbar << ": " << matWeigths(kbar,ktilde) << "\n";
}

std::complex<double> ThreadWrapper::gamma_oneD_spsp(double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar){
    // Watch out: do not use _Gk() because operator() involves static matrix in FunctorBuildGk and 
    // therefore generates data race when more then one thread used.
    std::complex<double> lower_level=0.0+0.0*im;
    for (int wttilde=0; wttilde<_Gk._size; wttilde++){
        for (size_t qttilde=0; qttilde<_Gk._kArr_l.size(); qttilde++){
            lower_level += buildGK(wtilde-_Gk._precomp_qn[wttilde],ktilde-_Gk._kArr_l[qttilde])[0]*buildGK(wbar-_Gk._precomp_qn[wttilde],kbar-_Gk._kArr_l[qttilde])[1];
            //lower_level += _Gk(wtilde-_Gk._precomp_qn[wttilde],ktilde-_Gk._kArr_l[qttilde])(0,0) * _Gk(wbar-_Gk._precomp_qn[wttilde],kbar-_Gk._kArr_l[qttilde])(1,1);
        }
    }
    
    lower_level *= SPINDEG*_Gk._u/(_Gk._beta*_Gk._Nk); /// Removed minus sign
    lower_level += 1.0;
    //cout << "gamma_oneD_spsp: " << _Gk._u/lower_level << endl;
    return _Gk._u/lower_level;
}

inline std::vector< std::complex<double> > ThreadWrapper::buildGK(std::complex<double> ik, double k){
    std::vector< std::complex<double> > GK = { 1.0/( ik + _Gk._mu - _Gk.epsk1D(k) - _Gk._u*_ndo_converged ), 1.0/( ik + _Gk._mu - _Gk.epsk1D(k) - _Gk._u*(1.0-_ndo_converged) ) }; // UP, DOWN
    return GK;
}

void ThreadWrapper::join_all(std::vector<std::thread>& grp){
    for (auto& thread : grp){
        if (thread.joinable()){
            thread.join();
        }
    }
}