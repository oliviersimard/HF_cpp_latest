#include "src/susceptibility_utils.h"
//#include "src/green_utils.h"
#include <iomanip>

using namespace std;

#define ONED

int main(int argc, char** argv){

    int Nomega=70;
    int Nk=50;
	double beta=5;
    double ndo_initial=0.6;
    int Niterations=100;
    double mu=0.;
    double ndo=ndo_initial;
    vector<complex<double> > Gup_k(Nomega);
    Susceptibility susObj;
    vector<double> kArr(Nk+1), kArr_l(Nk+1);
    for (int k=0; k<=Nk; k++){
        kArr[k] = -1.0*M_PI/2.0 + k*2.0*M_PI/(2.0*Nk);
        kArr_l[k] = -1.0*M_PI + k*2.0*M_PI/Nk;
    }

  
    for (double u=1.0; u<8; u+=0.2) {
        mu=u/2.;
        ndo=ndo_initial;

        Hubbard::FunctorBuildGk u_ndo_c(mu,beta,u,ndo,kArr,kArr_l,Niterations,Nk,Gup_k);
        #ifdef ONED
        complex<double> test;
        if (VERBOSE > 0) cout << "First 1D: " << u_ndo_c << endl;
        u_ndo_c.get_ndo_1D();
        cout << "After 1D: " << u_ndo_c << endl;
        test = susObj.chispsp_long_expr(u_ndo_c,Hubbard::K_1D(0.0,0.0+0.0*im));
        cout << "test here: " << test << endl;
        #else
        // //arma::Mat< complex<double> > test;
        if (VERBOSE > 0) cout << "First 2D: " << u_ndo_c << endl;
        u_ndo_c.get_ndo_2D();
        cout << "After 2D: " << u_ndo_c << endl;
        // //test = someFunction(u_ndo_c,kArr,Nomega);
        #endif    
    }

    return 0;
}
