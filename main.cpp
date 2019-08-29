#include "src/susceptibility_utils.h"
#include <iomanip>
#include <sys/stat.h>
//#include <omp.h>

inline bool file_exists(const std::string&);
double funct(double x){
    return x*x;
}

using namespace std;

// #define ONED

int main(int argc, char** argv){

    cout << "Number of threads: " << thread::hardware_concurrency() << '\n';

    const string filename("params.json");
    Json_utils JsonObj;
    MembCarrier params = JsonObj.JSONLoading(filename); // Loading file content into container.

    const int Nomega=*(params.int_ptr);
    const int Nk=*(params.int_ptr+3);
	const double beta_init=*(params.db_ptr+8),beta_step=*(params.db_ptr+7),beta_max=*(params.db_ptr+6);
    const double u_init=*(params.db_ptr+2),u_step=*(params.db_ptr+1),u_max=*(params.db_ptr);

    // Integrals integralsObj;
    // double result, result2;
    // result = integralsObj.integrate_simps(&funct,-2.0,2.0,0.0001);
    // result2 = integralsObj.integrate_trap(&funct,-2.0,2.0,0.0001);
    // cout << result << " " << result2 << endl;

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

    string testStr("_sum_for_bubb_real_part_minus_in_bubble"); // Should be "" when not testing. Adapt it otherwise (appends at end of every filenames.)
    string frontEnd("chiAM/"); // The folder in data/ containing the data.

    #ifdef ONED
    string fileOutput("data/"+frontEnd+"TvsU_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
    string fileOutputGtau("data/"+frontEnd+"Gtau_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
    string fileOutputMeanChi("data/"+frontEnd+"Mean_1D_chi_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
    string fileOutputMeanChi0("data/"+frontEnd+"Mean_1D_chio_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
    #else
    string fileOutput("data/"+frontEnd+"TvsU_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
    string fileOutputGtau("data/"+frontEnd+"Gtau_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
    string fileOutputMeanChi("data/"+frontEnd+"Mean_2D_chi_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
    string fileOutputMeanChi0("data/"+frontEnd+"Mean_2D_chio_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
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
    ofstream outputFileMeanChi;
    ofstream outputFileMeanChi0;
    for (double beta=beta_init; beta<=beta_max; beta+=beta_step){
        // Filenames for the file outputs of the susceptibilities chi and chi0.

        #ifdef ONED
        string fileOutputChi("data/"+frontEnd+"Chisp_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
        string fileOutputChi0("data/"+frontEnd+"Chi0_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
        #else
        string fileOutputChi("data/"+frontEnd+"Chisp_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
        string fileOutputChi0("data/"+frontEnd+"Chi0_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
        #endif

        for (double u=u_init; u<=u_max; u+=u_step) {
            mu=u/2.0;
            ndo=ndo_initial;
            vector<double> Gtau;
            Hubbard::FunctorBuildGk u_ndo_c(mu,beta,u,ndo,kArr,kArr_l,Niterations,Nk,Gup_k);
            #ifdef ONED
            if (VERBOSE > 0) cout << "First 1D: " << u_ndo_c << endl;
            /* Getting the selfconsistency part done. */
            u_ndo_c.get_ndo_1D();
            /* Printing result after selfconsistency. */
            cout << "After 1D: " << u_ndo_c << endl;
            // saveGF_grid(string("data/iwn_k_grid_GF"+testStr),u_ndo_c); // Saves G(tau) for each k-value. Eventually add as optional!!
            /* Computing ladder susceptibility diagram. */
            // complex<double> susLadder;
            // susLadder = susObj.chispsp_long_expr(u_ndo_c,Hubbard::K_1D(0.0,0.0+0.0*im));
            // cout << "Ladder susceptibility: " << susLadder << endl;
            /* Getting chi and chi0 printed into files (4). */
            tuple< complex<double>, complex<double> > chis;
            chis=susObj.get_chi_1D(u_ndo_c,fileOutputChi,fileOutputChi0); // Prints out chi and chi0 into files. Eventually add a boolean value to select this as option in params.json.
            outputFileMeanChi.open(fileOutputMeanChi, ofstream::out | ofstream::app);
            outputFileMeanChi0.open(fileOutputMeanChi0, ofstream::out | ofstream::app);
            outputFileMeanChi << get<0>(chis).real() << " ";
            outputFileMeanChi0 << get<1>(chis).real() << " ";
            outputFileMeanChi.close();
            outputFileMeanChi0.close();
            /* Getting G(\tau) for each value */
            Gtau = fftObj.get_gtau1D(u_ndo_c); // This saves the k-averaged G(tau) Green‘s function.
            outputFileGtau.open(fileOutputGtau, ofstream::out | ofstream::app);
            for (double el : Gtau){
                outputFileGtau << el << " ";
            }
            outputFileGtau << "\n";
            outputFileGtau.close();
            #else
            // complex<double> testSup;
            if (VERBOSE > 0) cout << "First 2D: " << u_ndo_c << endl;
            u_ndo_c.get_ndo_2D();
            cout << "After 2D: " << u_ndo_c << endl;

            /* Computing the 2D susceptibility */
            // testSup = susObj.chisp(u_ndo_c,Hubbard::K_2D(M_PI,M_PI,0.0+0.0*im));
            // cout << "test susceptibility: " << testSup << endl;
            tuple< complex<double>, complex<double> > chis;
            chis=susObj.get_chi_2D(u_ndo_c,fileOutputChi,fileOutputChi0);
            outputFileMeanChi.open(fileOutputMeanChi, ofstream::out | ofstream::app);
            outputFileMeanChi0.open(fileOutputMeanChi0, ofstream::out | ofstream::app);
            outputFileMeanChi << get<0>(chis) << " ";
            outputFileMeanChi0 << get<1>(chis) << " ";
            outputFileMeanChi.close();
            outputFileMeanChi0.close();

            /* Getting G(\tau) for each value */
            Gtau = fftObj.get_gtau2D(u_ndo_c);
            outputFileGtau.open(fileOutputGtau, ofstream::out | ofstream::app);
            for (double el : Gtau){
                outputFileGtau << el << " ";
            }
            outputFileGtau << "\n";
            outputFileGtau.close();
            #endif
            
            // outputFile.open(fileOutput, ofstream::out | ofstream::app);
            // outputFile << u_ndo_c.get_ndo() << " ";
            // outputFile.close();
        }
        // outputFile.open(fileOutput, ofstream::out | ofstream::app);
        // outputFile << "\n";
        // outputFile.close();
        outputFileMeanChi.open(fileOutputMeanChi, ofstream::out | ofstream::app);
        outputFileMeanChi << "\n";
        outputFileMeanChi.close();
        outputFileMeanChi0.open(fileOutputMeanChi0, ofstream::out | ofstream::app);
        outputFileMeanChi0 << "\n";
        outputFileMeanChi0.close();
    }

    return 0;
}


inline bool file_exists (const std::string& filename) {
  struct stat buffer;   
  return (stat(filename.c_str(), &buffer) == 0); 
}