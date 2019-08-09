#include "src/susceptibility_utils.h"
#include <iomanip>
#include <sys/stat.h>

inline bool file_exists(const std::string&);

using namespace std;

#define ONED

int main(int argc, char** argv){

    const string filename("params.json");
    Json_utils JsonObj;
    MembCarrier params = JsonObj.JSONLoading(filename); // Loading file content into container.

    const int Nomega=*(params.int_ptr);
    const int Nk=*(params.int_ptr+3);
	const double beta_init=*(params.db_ptr+8),beta_step=*(params.db_ptr+7),beta_max=*(params.db_ptr+6);
    const double u_init=*(params.db_ptr+2),u_step=*(params.db_ptr+1),u_max=*(params.db_ptr);

    const double ndo_initial=0.6;
    const int Niterations=*(params.int_ptr+1);
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
    string fileOutput("data/TvsU_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    string fileOutputGtau("data/Gtau_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    #else
    string fileOutput("data/TvsU_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    string fileOutputGtau("data/Gtau_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    string fileOutputChi("data/Chisp_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
    #endif
    // Testing if the files already exist.
    if (file_exists(fileOutput)){
        cout << "file: " << fileOutput << " already exists!" << endl;
        exit(0);
    }
    if (file_exists(fileOutputGtau)){
        cout << "file: " << fileOutputGtau << " already exists!" << endl;
        exit(0);
    }
    ofstream outputFile;
    ofstream outputFileGtau;
    ofstream outputFileChi;
    ofstream outputFileChi0;
    for (double beta=beta_init; beta<=beta_max; beta+=beta_step){

        string fileOutputChi("data/Chisp_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
        string fileOutputChi0("data/Chi0_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+".dat");
        for (double u=u_init; u<=u_max; u+=u_step) {
            mu=u/2.0;
            ndo=ndo_initial;
            std::vector<double> Gtau;
            Hubbard::FunctorBuildGk u_ndo_c(mu,beta,u,ndo,kArr,kArr_l,Niterations,Nk,Gup_k);
            #ifdef ONED
            complex<double> Suscep, Suscep0;
            if (VERBOSE > 0) cout << "First 1D: " << u_ndo_c << endl;
            /* Getting the selfconsistency part done. */
            u_ndo_c.get_ndo_1D();
            /* Printing result after selfconsistency. */
            cout << "After 1D: " << u_ndo_c << endl;
            saveGF_grid(string("data/iwn_k_grid_GF"),u_ndo_c);
            /* Computing ladder susceptibility diagram. */
            for (int k=0; k<=Nk; k++){
                Suscep = susObj.chisp(u_ndo_c,Hubbard::K_1D(kArr_l[k],0.0+0.0*im)); // Bosonic Matsubara frequency!
                Suscep0 = susObj.chi0(u_ndo_c,Hubbard::K_1D(kArr_l[k],0.0+0.0*im)); // Bosonic Matsubara frequency!
                /* Saving imaginary part of the susceptibility */
                outputFileChi.open(fileOutputChi, ofstream::out | ofstream::app);
                outputFileChi << 1.0*Suscep.imag() << " ";
                outputFileChi.close();
                // cout << "susceptibility at k: " << kArr_l[k] << " is " << Suscep << endl;
                outputFileChi0.open(fileOutputChi0, ofstream::out | ofstream::app);
                outputFileChi0 << Suscep0.imag() << " ";
                outputFileChi0.close();
            }
            cout << "susceptibility: " << Suscep << endl;
            outputFileChi.open(fileOutputChi, ofstream::out | ofstream::app);
            outputFileChi << "\n";
            outputFileChi.close();
            outputFileChi0.open(fileOutputChi0, ofstream::out | ofstream::app);
            outputFileChi0 << "\n";
            outputFileChi0.close();
    
            /* Getting G(\tau) for each value */
            Gtau = fftObj.get_gtau1D(u_ndo_c);
            outputFileGtau.open(fileOutputGtau, ofstream::out | ofstream::app);
            for (double el : Gtau){
                outputFileGtau << el << " ";
            }
            outputFileGtau << "\n";
            outputFileGtau.close();
            #else
            complex<double> testSup;
            if (VERBOSE > 0) cout << "First 2D: " << u_ndo_c << endl;
            u_ndo_c.get_ndo_2D();
            cout << "After 2D: " << u_ndo_c << endl;

            /* Computing the 2D susceptibility */
            testSup = susObj.chisp(u_ndo_c,Hubbard::K_2D(M_PI,M_PI,0.0+0.0*im));
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
    }

    return 0;
}


inline bool file_exists (const std::string& filename) {
  struct stat buffer;   
  return (stat(filename.c_str(), &buffer) == 0); 
}