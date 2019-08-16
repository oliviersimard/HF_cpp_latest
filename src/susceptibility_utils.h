#ifndef SUSCEPTIBILITY_H_
#define SUSCEPTIBILITY_H_

#include "fft.h"
#include "integral_utils.h"
#include "thread_utils.h"

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

};

#endif /* SUSCEPTIBILITY_H_ */