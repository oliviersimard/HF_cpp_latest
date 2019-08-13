#ifndef SUSCEPTIBILITY_H_
#define SUSCEPTIBILITY_H_

#include "fft.h"
#include "integral_utils.h"

class Susceptibility{ // The current-current susceptibility uses the +U in gamma, multiplied by the currents thereafter.
    public:
        std::complex<double> gamma_oneD_spsp(Hubbard::FunctorBuildGk&,Hubbard::K_1D,Hubbard::K_1D,Hubbard::K_1D) const;
        std::complex<double> gamma_oneD_spsp(Hubbard::FunctorBuildGk&,double,std::complex<double>,double,std::complex<double>,Hubbard::K_1D) const;
        std::complex<double> chispsp(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;   
        std::complex<double> chispsp_long_expr(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;
        std::complex<double> chisp(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;
        std::complex<double> chisp(Hubbard::FunctorBuildGk&,Hubbard::K_2D) const;
        std::complex<double> chi0(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;
        std::complex<double> chi0(Hubbard::FunctorBuildGk&,Hubbard::K_2D) const;
        std::complex<double> gamma_oneD_spsp_full_lower(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D ktilde,Hubbard::K_1D kbar,Hubbard::K_1D qtilde,Hubbard::K_1D q) const;
        std::complex<double> gamma_oneD_spsp_full_middle(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D kbar,Hubbard::K_1D q) const;
        std::complex<double> chisp_full(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const;
        void get_chi_1D(Hubbard::FunctorBuildGk& Gk, std::string filename_chi, std::string filename_chi0);

        // inline void wait(int seconds){
        //     boost::this_thread::sleep_for(boost::chrono::seconds(seconds));
        // }

};

// template <typename Lockable>
// class strict_lock  { // For external locking when parallelizing code! Not used for now.
// public:
//     typedef Lockable lockable_type;

//     explicit strict_lock(lockable_type& obj) : obj_(obj) {
//         obj.lock(); // locks on construction
//     }
//     strict_lock() = delete;
//     strict_lock(const strict_lock&) = delete;
//     strict_lock& operator=(const strict_lock&) = delete;

//     ~strict_lock() { obj_.unlock(); } //  unlocks on destruction 

//     bool owns_lock(const lockable_type* l) const noexcept // strict lockers specific function 
//     {
//       return l == &obj_;
//     }
// private:
//     lockable_type& obj_;
// };

// class ThreadFunctor{
//     public:
//         ThreadFunctor(size_t ktilde,size_t kbar,std::complex<double>& upper_level,std::complex<double>& tmp_val_kbar,Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) : 
//         _ktilde(ktilde), _kbar(kbar), _upper_level(upper_level), _tmp_val_kbar(tmp_val_kbar), _Gk(Gk), _q(q) {
//             Susceptibility susObj;
//             this->_susObj=susObj;
//         }
//         ~ThreadFunctor(){};
//         void operator()(){
//             for (int wtilde=0; wtilde<_Gk._size; wtilde++){
//                 for (int wbar=0; wbar<_Gk._size; wbar++){
//                     std::complex<double> tmp_val = _Gk(
//                             _Gk._precomp_wn[wtilde],_Gk._kArr_l[_ktilde]
//                             )(0,0)*_Gk(
//                             _Gk._precomp_wn[wtilde]-_q._iwn,_Gk._kArr_l[_ktilde]-_q._qx
//                             )(0,0)*_susObj.gamma_oneD_spsp( _Gk, _Gk._kArr_l[_ktilde], _Gk._precomp_wn[wtilde], _Gk._kArr_l[_kbar], _Gk._precomp_wn[wbar], _q )*_Gk(
//                             _Gk._precomp_wn[wbar]+_q._iwn,_Gk._kArr_l[_kbar]+_q._qx
//                             )(1,1)*_Gk(
//                             _Gk._precomp_wn[wbar],_Gk._kArr_l[_kbar]
//                             )(1,1);
//                     boost::lock_guard<boost::mutex> guard(_mtx); // Internally locked
//                     _upper_level += tmp_val;
//                     _tmp_val_kbar += tmp_val;
//                 } 
//             }
//         }
//         void lock() {
//             _mtx.lock();
//         }
//         void unlock() {
//             _mtx.unlock();
//         }
//     private:
//         boost::mutex _mtx;
//         size_t _ktilde, _kbar;
//         std::complex<double> _upper_level, _tmp_val_kbar;
//         Hubbard::FunctorBuildGk _Gk;
//         Hubbard::K_1D _q;
//         Susceptibility _susObj;
// };

#endif /* SUSCEPTIBILITY_H_ */