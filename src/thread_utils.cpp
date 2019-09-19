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
    lower_level *= 1.0*_Gk._u/(_Gk._beta*_Gk._Nk); /// Removed minus sign
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

ThreadWrapper::ThreadWrapper(Hubbard::FunctorBuildGk Gk,Hubbard::K_2D& q,double ndo_converged) : _q(q){
    this->_ndo_converged=ndo_converged;
    this->_Gk=Gk;
}

void ThreadWrapper::operator()(int ktilde, int kbar, double beta){
    std::complex<double> tmp_val_kt_kb(0.0,0.0);
    std::complex<double> tmp_val_weights(0.0,0.0);
    std::complex<double> tmp_val_tot_sus(0.0,0.0);
    // cout << beta << " " << _Gk._kArr_l[ktilde] << " " << _q._iwn << " " << _q._qx << " " << _Gk._kArr_l[kbar] << endl;
    for (int wtilde=0; wtilde<_Gk._size; wtilde++){
        for (int wbar=0; wbar<_Gk._size; wbar++){
            std::complex<double> tmp_val_kt_kb_tmp = gamma_oneD_spsp(_Gk._kArr_l[ktilde],_Gk._precomp_wn[wtilde],_Gk._kArr_l[kbar],_Gk._precomp_wn[wbar]);
            std::complex<double> tmp_val_weights_tmp = buildGK1D(
                                    _Gk._precomp_wn[wtilde],_Gk._kArr_l[ktilde]
                                    )[0]*buildGK1D(
                                    _Gk._precomp_wn[wtilde]-_q._iwn,_Gk._kArr_l[ktilde]-_q._qx
                                    )[0]*buildGK1D(
                                    _Gk._precomp_wn[wbar]-_q._iwn,_Gk._kArr_l[kbar]-_q._qx
                                    )[1]*buildGK1D(
                                    _Gk._precomp_wn[wbar],_Gk._kArr_l[kbar]
                                    )[1]; 
            tmp_val_weights += tmp_val_weights_tmp;
            tmp_val_kt_kb += tmp_val_kt_kb_tmp;
            tmp_val_tot_sus += tmp_val_kt_kb_tmp*tmp_val_weights_tmp;
            if ((wtilde==0) && (wbar==0)){
                std::cout << "Thread id: " << std::this_thread::get_id() << std::endl;
                std::cout << "ktilde: " << _Gk._kArr_l[ktilde] << std::endl;
                std::cout << "kbar: " << _Gk._kArr_l[kbar] << std::endl;
                std::cout << "gamma_oneD_spsp: " << tmp_val_kt_kb << std::endl;
                std::cout << "weights: " << tmp_val_weights << std::endl;
            }
        } 
    }
    // lock_guard<mutex> guard(mutx); 
    matGamma(kbar,ktilde) = tmp_val_kt_kb; // These matrices are static variables.
    matWeigths(kbar,ktilde) = tmp_val_weights;
    matTotSus(kbar,ktilde) = 1.0/(_Gk._beta)/(_Gk._beta)*tmp_val_tot_sus; // This gives the total susceptibility resolved in k-space. Summation performed on beta only.
    //cout << "Gamma for " << "ktilde " << ktilde << " and kbar " << kbar << ": " << matGamma(kbar,ktilde) << "\n";
    //cout << "Weigths for " << "ktilde " << ktilde << " and kbar " << kbar << ": " << matWeigths(kbar,ktilde) << "\n";
}

void ThreadWrapper::operator()(int kbarx_m_tildex, int kbary_m_tildey){ // no beta to overload operator() function.
/* In 2D, calculating the weights implies computing whilst dismissing the translational invariance of the vertex function (Gamma). */
    std::complex<double> tmp_val_kt_kb(0.0,0.0);
    // cout << beta << " " << _Gk._kArr_l[ktilde] << " " << _q._iwn << " " << _q._qx << " " << _Gk._kArr_l[kbar] << endl;
    for (int wtilde=0; wtilde<_Gk._size; wtilde++){
        for (int wbar=0; wbar<_Gk._size; wbar++){
            tmp_val_kt_kb += gamma_twoD_spsp(_Gk._kArr_l[kbarx_m_tildex],_Gk._kArr_l[kbary_m_tildey],_Gk._precomp_wn[wtilde],_Gk._precomp_wn[wbar]);
            if ((wtilde==0) && (wbar==0)){
                std::cout << "Thread id: " << std::this_thread::get_id() << std::endl;
                std::cout << "kbarx_m_tildex: " << _Gk._kArr_l[kbarx_m_tildex] << std::endl;
                std::cout << "kbary_m_tildey: " << _Gk._kArr_l[kbary_m_tildey] << std::endl;
                std::cout << "gamma_oneD_spsp: " << tmp_val_kt_kb << std::endl;
            }
        } 
    }
    // lock_guard<mutex> guard(mutx);
    matGamma(kbary_m_tildey,kbarx_m_tildex) = tmp_val_kt_kb; // These matrices are static variables.
}

std::complex<double> ThreadWrapper::gamma_oneD_spsp(double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar){
    // Watch out: do not use _Gk() because operator() involves static matrix in FunctorbuildGk and 
    // therefore generates data race when more then one thread used.
    std::complex<double> lower_level=0.0+0.0*im;
    for (int wttilde=0; wttilde<_Gk._size; wttilde++){
        for (size_t qttilde=0; qttilde<_Gk._kArr_l.size(); qttilde++){
            lower_level += buildGK1D(wtilde-_Gk._precomp_qn[wttilde],ktilde-_Gk._kArr_l[qttilde])[0]*buildGK1D(wbar-_Gk._precomp_qn[wttilde],kbar-_Gk._kArr_l[qttilde])[1];
            //lower_level += _Gk(wtilde-_Gk._precomp_qn[wttilde],ktilde-_Gk._kArr_l[qttilde])(0,0) * _Gk(wbar-_Gk._precomp_qn[wttilde],kbar-_Gk._kArr_l[qttilde])(1,1);
        }
    }
    
    lower_level *= SPINDEG*_Gk._u/(_Gk._beta*_Gk._Nk); /// Removed minus sign
    lower_level += 1.0;
    //cout << "gamma_oneD_spsp: " << _Gk._u/lower_level << endl;
    return _Gk._u/lower_level;
}

std::complex<double> ThreadWrapper::gamma_twoD_spsp(double kbarx_m_tildex,double kbary_m_tildey,std::complex<double> wtilde,std::complex<double> wbar){
    std::complex<double> lower_level(0.0,0.0);
    for (int wttilde=0; wttilde<_Gk._size; wttilde++){
        for (size_t qttildey=0; qttildey<_Gk._kArr_l.size(); qttildey++){
            for (size_t qttildex=0; qttildex<_Gk._kArr_l.size(); qttildex++){ // the change of variable only applies to k-space, due to periodicity modulo 2pi.
                lower_level += buildGK2D(-(wtilde-_Gk._precomp_qn[wttilde]),-_Gk._kArr_l[qttildex],-_Gk._kArr_l[qttildey])[0]*buildGK2D(-(wbar-_Gk._precomp_qn[wttilde]),-(_Gk._kArr_l[qttildex]+kbarx_m_tildex),-(_Gk._kArr_l[qttildey]+kbary_m_tildey))[1];
            }
        }
    }
    lower_level *= SPINDEG*_Gk._u/(_Gk._beta*_Gk._Nk*_Gk._Nk); //factor 2 for the spin and minus sign added
    lower_level += 1.0;

    return _Gk._u/lower_level;
}

inline std::vector< std::complex<double> > ThreadWrapper::buildGK1D(std::complex<double> ik, double k){
    std::vector< std::complex<double> > GK = { 1.0/( ik + _Gk._mu - _Gk.epsk1D(k) - _Gk._u*_ndo_converged ), 1.0/( ik + _Gk._mu - _Gk.epsk1D(k) - _Gk._u*(1.0-_ndo_converged) ) }; // UP, DOWN
    return GK;
}

inline std::vector< std::complex<double> > ThreadWrapper::buildGK2D(std::complex<double> ik, double kx, double ky){
    std::vector< std::complex<double> > GK = { 1.0/( ik + _Gk._mu - _Gk.epsk2D(kx,ky) - _Gk._u*_ndo_converged ), 1.0/( ik + _Gk._mu - _Gk.epsk2D(kx,ky) - _Gk._u*(1.0-_ndo_converged) ) }; // UP, DOWN
    return GK;
}

void ThreadWrapper::join_all(std::vector<std::thread>& grp){
    for (auto& thread : grp){
        if (thread.joinable()){
            thread.join();
        }
    }
}