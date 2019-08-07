#include "src/susceptibility_utils.h"
#include <iomanip>

using namespace std;

#define ONED

int main(int argc, char** argv){

    int Nomega=400;
    int Nk=200;
	double beta_init=5.0,beta_step=5.0,beta_max=100.0;
    double u_init=0.0,u_step=0.1,u_max=8.0;
    double ndo_initial=0.6;
    int Niterations=100;
    double mu=0.0;
    double ndo;
    vector< complex<double> > Gup_k(Nomega,0.0);
    Susceptibility susObj;
    FFT fftObj;
    vector<double> kArr(Nk+1), kArr_l(Nk+1);
    for (int k=0; k<=Nk; k++){
        kArr[k] = -1.0*M_PI/2.0 + k*2.0*M_PI/(2.0*Nk);
        kArr_l[k] = -1.0*M_PI + k*2.0*M_PI/Nk;
    }

    #ifdef ONED
    string fileOutput("TvsU_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    string fileOutputGtau("Gtau_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    string fileOutputChi("Chispsp_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    #else
    string fileOutput("TvsU_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    string fileOutputGtau("Gtau_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    string fileOutputChi("Chisp_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    #endif
    ofstream outputFile;
    ofstream outputFileGtau;
    ofstream outputFileChi;
    for (double beta=beta_init; beta<=beta_max; beta+=beta_step){
        
        for (double u=u_init; u<=u_max; u+=u_step) {
            mu=u/2.;
            ndo=ndo_initial;
            std::vector<double> Gtau;
            Hubbard::FunctorBuildGk u_ndo_c(mu,beta,u,ndo,kArr,kArr_l,Niterations,Nk,Gup_k);
            #ifdef ONED
            complex<double> testSup;
            if (VERBOSE > 0) cout << "First 1D: " << u_ndo_c << endl;
            /* Getting the selfconsistency part done. */
            u_ndo_c.get_ndo_1D();
            /* Printing result after selfconsistency. */
            cout << "After 1D: " << u_ndo_c << endl;
            /* Computing ladder susceptibility diagram. */
            //testSup = susObj.chisp(u_ndo_c,Hubbard::K_1D(0.0,0.0+0.0*im));
            testSup = susObj.chisp(u_ndo_c,Hubbard::K_1D(0.0,0.0+0.0*im));
            cout << "test susceptibility: " << testSup << endl;

            /* Getting G(\tau) for each value */
            Gtau = fftObj.get_gtau1D(u_ndo_c);
            outputFileGtau.open(fileOutputGtau, ofstream::out | ofstream::app);
            for (double el : Gtau){
                outputFileGtau << el << " ";
            }
            outputFileGtau << "\n";
            outputFileGtau.close();

            /* Saving imaginary part of the susceptibility */
            outputFileChi.open(fileOutputChi, ofstream::out | ofstream::app);
            outputFileChi << testSup.imag() << " ";
            outputFileChi.close();
            #else
            complex<double> testSup;
            if (VERBOSE > 0) cout << "First 2D: " << u_ndo_c << endl;
            u_ndo_c.get_ndo_2D();
            cout << "After 2D: " << u_ndo_c << endl;

            /* Computing the 2D susceptibility */
            testSup = susObj.chisp(u_ndo_c,Hubbard::K_2D(0.0,0.0,0.0+0.0*im));
            cout << "test susceptibility: " << testSup << endl;

            /* Getting G(\tau) for each value */
            Gtau = fftObj.get_gtau2D(u_ndo_c);
            outputFileGtau.open(fileOutputGtau, ofstream::out | ofstream::app);
            for (double el : Gtau){
                outputFileGtau << el << " ";
            }
            outputFileGtau << "\n";
            outputFileGtau.close();

            /* Saving imaginary part of the susceptibility */
            outputFileChi.open(fileOutputChi, ofstream::out | ofstream::app);
            outputFileChi << testSup.imag() << " ";
            outputFileChi.close();
            #endif
            
            outputFile.open(fileOutput, ofstream::out | ofstream::app);
            outputFile << u_ndo_c.get_ndo() << " ";
            outputFile.close();
        }
        outputFile.open(fileOutput, ofstream::out | ofstream::app);
        outputFile << "\n";
        outputFile.close();

        outputFileChi.open(fileOutputChi, ofstream::out | ofstream::app);
        outputFileChi << "\n";
        outputFileChi.close();
    }

    return 0;
}
