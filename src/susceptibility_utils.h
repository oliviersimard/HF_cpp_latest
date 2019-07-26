#ifndef SUSCEPTIBILITY_H_
#define SUSCEPTIBILITY_H_

#include "green_utils.h"
#include <string>

class Susceptibility{ // The current-current susceptibility uses the +U in gamma, multiplied by the currents thereafter.
    public:
        std::complex<double> gamma_oneD_spsp(Hubbard::FunctorBuildGk&,Hubbard::K_1D,Hubbard::K_1D,Hubbard::K_1D) const;
        std::complex<double> gamma_oneD_spsp(Hubbard::FunctorBuildGk&,double,std::complex<double>,double,std::complex<double>,Hubbard::K_1D) const;
        std::complex<double> chispsp(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;   
        std::complex<double> chispsp_long_expr(Hubbard::FunctorBuildGk&,Hubbard::K_1D) const;
};

#endif /* SUSCEPTIBILITY_H_ */