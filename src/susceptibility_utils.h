#ifndef SUSCEPTIBILITY_H_
#define SUSCEPTIBILITY_H_


#include "fft.h"
#include "integral_utils.h"
#include "thread_utils.h"
#include <tuple>
#include <algorithm> // For std::move

class Susceptibility{ // The current-current susceptibility uses the +U in gamma, multiplied by the currents thereafter.
    public:
        std::complex<double> gamma_oneD_spsp(Hubbard::FunctorBuildGk&,Hubbard::K_1D,Hubbard::K_1D,Hubbard::K_1D) const;
        std::complex<double> gamma_oneD_spsp(Hubbard::FunctorBuildGk&,double,std::complex<double>,double,std::complex<double>) const;
        std::complex<double> gamma_twoD_spsp(Hubbard::FunctorBuildGk& Gk,double kbarx_m_tildex,double kbary_m_tildey,std::complex<double> wtilde,std::complex<double> wbar) const;
        std::complex<double> chispsp(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;   
        std::complex<double> chispsp_long_expr(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;
        std::complex<double> chisp(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;
        std::complex<double> chisp(Hubbard::FunctorBuildGk&,Hubbard::K_2D) const;
        std::complex<double> chi0(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;
        std::complex<double> chi0(Hubbard::FunctorBuildGk&,Hubbard::K_2D) const;
        std::complex<double> gamma_oneD_spsp_full_lower(Hubbard::FunctorBuildGk& Gk,double kp,double kbar,std::complex<double> iknp,std::complex<double> wbar) const;
        std::complex<double> gamma_twoD_spsp_full_lower(Hubbard::FunctorBuildGk& Gk,double kpx,double kpy,double kbarx,double kbary,std::complex<double> iknp,std::complex<double> wbar) const;
        std::tuple< std::complex<double>, std::complex<double>, std::complex<double> > gamma_oneD_spsp_full_middle_plotting(Hubbard::FunctorBuildGk& Gk,double kbar,double ktilde,std::complex<double> wbar,std::complex<double> wtilde,Hubbard::K_1D q) const;
        std::tuple< std::complex<double>,std::complex<double>,std::complex<double> > gamma_twoD_spsp_full_middle_plotting(Hubbard::FunctorBuildGk& Gk,double kbarx_m_tildex,double kbary_m_tildey,std::complex<double> wbar,std::complex<double> wtilde,Hubbard::K_2D q) const;
        std::tuple< std::complex<double>,std::complex<double>,std::complex<double> > gamma_oneD_jj_full_middle_plotting(Hubbard::FunctorBuildGk& Gk,double kbar,double ktilde,std::complex<double> wbar,std::complex<double> wtilde,Hubbard::K_1D q) const;
        std::complex<double> chisp_full(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const;
        std::tuple< std::complex<double>, std::complex<double> > get_chi_1D(Hubbard::FunctorBuildGk& Gk, std::string filename_chi, std::string filename_chi0);
        std::tuple< std::complex<double>, std::complex<double> > get_chi_2D(Hubbard::FunctorBuildGk& Gk, std::string filename_chi, std::string filename_chi0);
        std::tuple< std::complex<double>, std::complex<double> > gamma_oneD_spsp_plotting(Hubbard::FunctorBuildGk& Gk,double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar,Hubbard::K_1D q) const;
        std::tuple< std::complex<double>, std::complex<double> > gamma_oneD_jj_plotting(Hubbard::FunctorBuildGk& Gk,double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar,Hubbard::K_1D q) const;
        std::tuple< std::complex<double>, std::complex<double> > gamma_twoD_spsp_plotting(Hubbard::FunctorBuildGk& Gk,double ktildex_m_barx,double ktildey_m_bary,std::complex<double> wtilde,std::complex<double> wbar) const;
        std::tuple< std::complex<double>, std::complex<double> > gamma_twoD_jj_plotting(Hubbard::FunctorBuildGk& Gk,double kbarx_m_tildex,double kbary_m_tildey,std::complex<double> wtilde,std::complex<double> wbar) const;
        std::tuple< std::complex<double>, std::complex<double> > gamma_oneD_spsp_crossed_plotting(Hubbard::FunctorBuildGk& Gk,double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar,Hubbard::K_1D q) const;
};

#endif /* SUSCEPTIBILITY_H_ */