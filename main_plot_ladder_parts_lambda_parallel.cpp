#include "src/susceptibility_utils.h"
#include <iomanip>
#include <sys/stat.h>
//#include <omp.h>

using namespace std;

inline bool file_exists(const string&);
void getWeights(Hubbard::FunctorBuildGk& Gk, Hubbard::K_2D qq, vector<double>& kArr_l,double beta,double Nomega,ofstream& fileObj,string& filename);
void join_all(vector<thread>&);

#define ONED

#ifdef PARALLEL
/* Remember that armadillo is column-major. Useful for parallel version. */
static arma::Mat< complex<double> > matGamma;
static arma::Mat< complex<double> > matWeigths;

// class ThreadWrapper{
//     public:
//         ThreadWrapper(Hubbard::FunctorBuildGk& Gk, Hubbard::K_1D& q, arma::cx_dmat::iterator matPtr, arma::cx_dmat::iterator matWPtr);
//         ThreadWrapper(Hubbard::FunctorBuildGk& Gk, Hubbard::K_1D& q);
//         ~ThreadWrapper(){};
//         void operator()(int ktilde, int kbar, double b);
//         complex<double> gamma_oneD_spsp(double ktilde,complex<double> wtilde,double kbar,complex<double> wbar);
//         void join_all(vector<thread>& grp);
//     private:
//         Hubbard::FunctorBuildGk& _Gk;
//         Hubbard::K_1D& _q;
//         arma::cx_dmat::iterator _ktb;
//         arma::cx_dmat::iterator _ktbW;
// };
#endif


int main(int argc, char** argv){

    cout << "Number of threads: " << thread::hardware_concurrency() << '\n';

    // Extracting some parameters from the JSON file.
    const string filename("params.json");
    Json_utils JsonObj;
    MembCarrier params = JsonObj.JSONLoading(filename); // Loading file content into container.

    const int Nomega=*(params.int_ptr);
    const int Nk=*(params.int_ptr+3);
	const double beta_init=*(params.db_ptr+8),beta_step=*(params.db_ptr+7),beta_max=*(params.db_ptr+6);
    const double u_init=*(params.db_ptr+2),u_step=*(params.db_ptr+1),u_max=*(params.db_ptr);
    const bool is_full=*(params.boo_ptr);
    //

    const double ndo_initial=0.6;
    const int Niterations=*(params.int_ptr+1);
    double mu=0.0;
    double ndo;
    vector< complex<double> > Gup_k(Nomega,0.0);
    vector<double> kArr(Nk+1), kArr_l(Nk+1);
    for (int k=0; k<=Nk; k++){
        kArr[k] = -1.0*M_PI/2.0 + k*2.0*M_PI/(2.0*Nk);
        kArr_l[k] = -1.0*M_PI + k*2.0*M_PI/Nk;
    }

    string testStr("_parallelized_first_fermionic_freq_minus_lower_bubble_spin_"+to_string(static_cast<int>(SPINDEG))+""); // Should be "" when not testing. Adapt it otherwise (appends at end of every filenames.)
    string frontEnd(""); // The folder in data/ containing the data.

    #ifdef ONED
    string fileOutput("data/"+frontEnd+"TvsU_1D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
    #else
    string fileOutput("data/"+frontEnd+"TvsU_2D_mapping_U_"+to_string(u_init)+"_"+to_string(u_step)+"_"+to_string(u_max)+"_beta_"+to_string(beta_init)+"_"+to_string(beta_step)+"_"+to_string(beta_max)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
    #endif
    // Testing if the files already exist.
    if (file_exists(fileOutput)){
        cout << "file: " << fileOutput << " already exists!" << endl;
        exit(0);
    }
    ofstream outputFile;
    ofstream outputFileChispspWeights;
    ofstream outputFileChispspGamma;
    ofstream outputFileChispspGammaBubble;
    for (double beta=beta_init; beta<=beta_max; beta+=beta_step){
        
        for (double u=u_init; u<=u_max; u+=u_step) {

            // Filenames for the file outputs of the susceptibilities chi and chi0.
            #ifdef ONED
            string fileOutputChispspWeigths;
            string fileOutputChispspGamma;
            string fileOutputChispspGammaBubble;
            if (!is_full){
                fileOutputChispspWeigths = "data/"+frontEnd+"ChispspWeights_1D_U_"+to_string(u)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat";
                fileOutputChispspGamma = "data/"+frontEnd+"ChispspGamma_1D_U_"+to_string(u)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat";
                fileOutputChispspGammaBubble = "data/"+frontEnd+"ChispspBubble_1D_U_"+to_string(u)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat";
            }
            else{
                fileOutputChispspWeigths = "data/"+frontEnd+"ChispspWeights_1D_U_"+to_string(u)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+"_full_.dat";
                fileOutputChispspGamma = "data/"+frontEnd+"ChispspGamma_1D_U_"+to_string(u)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+"_full_.dat";
                fileOutputChispspGammaBubble = "data/"+frontEnd+"ChispspBubble_1D_U_"+to_string(u)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+"_full_.dat";
            }
            #else
            string fileOutputChispspWeigths("data/"+frontEnd+"ChispspWeights_2D_U_"+to_string(u)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
            string fileOutputChispspGamma("data/"+frontEnd+"ChispspGamma_2D_U_"+to_string(u)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
            string fileOutputChispspGammaBubble("data/"+frontEnd+"ChispspBubble_2D_U_"+to_string(u)+"_beta_"+to_string(beta)+"_Nomega"+to_string(Nomega)+"_Nk"+to_string(Nk)+testStr+".dat");
            #endif

            mu=u/2.0;
            ndo=ndo_initial;
            vector<double> Gtau;
            Hubbard::FunctorBuildGk u_ndo_c(mu,beta,u,ndo,kArr,kArr_l,Niterations,Nk,Gup_k);
            #ifdef ONED
            if (VERBOSE > 0) cout << "First 1D: " << u_ndo_c << endl;
            /* Getting the selfconsistency part done. */
            u_ndo_c.get_ndo_1D();
            double ndo_converged = u_ndo_c.get_ndo();
            /* Printing result after selfconsistency. */
            cout << "After 1D: " << u_ndo_c << endl;
            /* Computing ladder susceptibility diagram. */
            Hubbard::K_1D qq(0.0,0.0+0.0*im); // photon 4-vector
            #ifdef PARALLEL
            unsigned int it = 0;
            unsigned int totsize = kArr_l.size()*kArr_l.size(); // Nk+1 * Nk+1
            matGamma = arma::cx_dmat(kArr_l.size(),kArr_l.size(),arma::fill::zeros);
            matWeigths = arma::cx_dmat(kArr_l.size(),kArr_l.size(),arma::fill::zeros);
            while (it<=totsize){
                if (totsize % NUM_THREADS == 0){
                    vector<thread> tt(NUM_THREADS);
                    for (unsigned int k=0; k<NUM_THREADS; k++){
                        unsigned int tmpit = it+k;
                        unsigned int lkt = static_cast<int>(floor(tmpit/kArr_l.size())); // Samples the rows
                        unsigned int lkb = (tmpit % kArr_l.size()); // Samples the columns
                        cout << "up lkt: " << lkt << " lkb: " << lkb << "\n";
                        thread t1([kArr_l,qq,lkt,lkb,Nomega,beta,u,Nk,ndo_converged](){ // u_ndo_c shouldn't be passed by reference; data race.
                            complex<double> tmp_val_kt_kb(0.0,0.0);
                            complex<double> tmp_val_weigths(0.0,0.0);
                            for (int wtilde=0; wtilde<Nomega; wtilde++){
                                for (int wbar=0; wbar<Nomega; wbar++){
                                    auto gamma_lambda_expr = [&](){
                                        complex<double> lower_level=0.0+0.0*im;
                                        for (int wttilde=0; wttilde<Nomega; wttilde++){
                                            for (size_t qttilde=0; qttilde<kArr_l.size(); qttilde++){
                                                //lower_level += u_ndo_c(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)-complex<double>(0.0,(2.0*wttilde)*M_PI/beta),kArr_l[lkt]-kArr_l[qttilde])(0,0)*u_ndo_c(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)-complex<double>(0.0,(2.0*wttilde)*M_PI/beta),kArr_l[lkb]-kArr_l[qttilde])(1,1);
                                                lower_level += (1.0/(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)-complex<double>(0.0,(2.0*wttilde)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkt]-kArr_l[qttilde])-u*ndo_converged)) * (1.0/(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)-complex<double>(0.0,(2.0*wttilde)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkb]-kArr_l[qttilde])-u*(1.0-ndo_converged)));
                                            }
                                        }
                                        lower_level *= -1.0*SPINDEG*u/(beta*Nk); /// Should it be devided by Nk or Nk+1?
                                        //lower_level += 1.0;
                                        //cout << "gamma_oneD_spsp: " << _Gk._u << " " << _Gk._Nk << endl;
                                        return lower_level;//_Gk._u/lower_level
                                    };
                                    auto gamma_lambda_val = gamma_lambda_expr();
                                    complex<double> weights_val = (1.0/(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkt])-u*ndo_converged)
                                        )*(1.0/(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)+qq._iwn+u/2.0+2.0*cos(kArr_l[lkt]+qq._qx)-u*ndo_converged)
                                        )*(1.0/(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)+qq._iwn+u/2.0+2.0*cos(kArr_l[lkb]+qq._qx)-u*(1.0-ndo_converged))
                                        )*(1.0/(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkb])-u*(1.0-ndo_converged))
                                        );
                                    //lock_guard<mutex> guard(mutx);
                                    tmp_val_kt_kb += gamma_lambda_val;
                                    tmp_val_weigths += weights_val;
                                    if ((wtilde==0) && (wbar==0)){
                                        cout << "Thread id: " << this_thread::get_id() << endl;
                                        cout << "ktilde: " << kArr_l[lkt] << endl;
                                        cout << "kbar: " << kArr_l[lkb] << endl;
                                        cout << "gamma_oneD_spsp: " << tmp_val_kt_kb << endl;
                                        cout << "weights: " << tmp_val_weigths << endl;
                                    }
                                } 
                            }
                            cout << "thread id operating: " << this_thread::get_id() << endl;
                        }); // END OF THREAD t1!!
                        tt[k]=move(t1);
                    }
                    join_all(tt);
                }
                else if (totsize % NUM_THREADS != 0){
                    if ( (totsize-it)<NUM_THREADS ){
                        unsigned int superk=totsize-it;
                        vector<thread> tt(superk);
                        for (unsigned int k=0; k<superk; k++){
                            unsigned int tmpit=it+k;
                            unsigned int lkt = static_cast<int>(floor(tmpit/kArr_l.size())); // Samples the rows
                            unsigned int lkb = (tmpit % kArr_l.size()); // Samples the columns
                            thread t1([kArr_l,qq,lkt,lkb,Nomega,beta,Nk,u,ndo_converged](){ // Shouldn
                                complex<double> tmp_val_kt_kb(0.0,0.0);
                                complex<double> tmp_val_weigths(0.0,0.0);
                                for (int wtilde=0; wtilde<Nomega; wtilde++){
                                    for (int wbar=0; wbar<Nomega; wbar++){
                                        auto gamma_lambda_expr = [&](){
                                            complex<double> lower_level=0.0+0.0*im;
                                            for (int wttilde=0; wttilde<Nomega; wttilde++){
                                                for (size_t qttilde=0; qttilde<kArr_l.size(); qttilde++){
                                                    lower_level += (1.0/(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)-complex<double>(0.0,(2.0*wttilde)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkt]-kArr_l[qttilde])-u*ndo_converged)) * (1.0/(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)-complex<double>(0.0,(2.0*wttilde)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkb]-kArr_l[qttilde])-u*(1.0-ndo_converged)));
                                                }
                                            }
                                            lower_level *= -1.0*SPINDEG*u/(beta*Nk); /// Removed minus sign
                                            //lower_level += 1.0;
                                            //cout << "gamma_oneD_spsp: " << _Gk._u << " " << _Gk._Nk << endl;
                                            return lower_level;//_Gk._u/lower_level
                                        };
                                        //lock_guard<mutex> guard(mutx);
                                        tmp_val_kt_kb += gamma_lambda_expr();
                                        tmp_val_weigths += (1.0/(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkt])-u*ndo_converged)
                                        )*(1.0/(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)+qq._iwn+u/2.0+2.0*cos(kArr_l[lkt]+qq._qx)-u*ndo_converged)
                                        )*(1.0/(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)+qq._iwn+u/2.0+2.0*cos(kArr_l[lkb]+qq._qx)-u*(1.0-ndo_converged))
                                        )*(1.0/(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkb])-u*(1.0-ndo_converged))
                                        );
                                        if ((wtilde==0) && (wbar==0)){
                                            cout << "Thread id: " << this_thread::get_id() << endl;
                                            cout << "ktilde: " << kArr_l[lkt] << endl;
                                            cout << "kbar: " << kArr_l[lkb] << endl;
                                            cout << "gamma_oneD_spsp: " << tmp_val_kt_kb << endl;
                                            cout << "weights: " << tmp_val_weigths << endl;
                                        }
                                    } 
                                }
                                cout << "thread id operating: " << this_thread::get_id() << endl;
                            }); // END OF THREAD t1!!
                            tt[k]=move(t1);

                        }
                        join_all(tt);
                    }
                    else{
                        vector<thread> tt(NUM_THREADS);
                        for (unsigned int k=0; k<NUM_THREADS; k++){
                            unsigned int tmpit=it+k;
                            unsigned int lkt = static_cast<int>(floor(tmpit/kArr_l.size())); // Samples the rows
                            unsigned int lkb = (tmpit % kArr_l.size()); // Samples the columns
                            cout << "down " << "lkt: " << lkt << " lkb: " << lkb << "\n";
                            thread t1([kArr_l,qq,lkt,lkb,Nomega,beta,Nk,u,ndo_converged](){
                                complex<double> tmp_val_kt_kb(0.0,0.0);
                                complex<double> tmp_val_weigths(0.0,0.0);
                                for (int wtilde=0; wtilde<Nomega; wtilde++){
                                    for (int wbar=0; wbar<Nomega; wbar++){
                                        auto gamma_lambda_expr = [&](){
                                            complex<double> lower_level=0.0+0.0*im;
                                            for (int wttilde=0; wttilde<Nomega; wttilde++){
                                                for (size_t qttilde=0; qttilde<kArr_l.size(); qttilde++){
                                                    lower_level += (1.0/(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)-complex<double>(0.0,(2.0*wttilde)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkt]-kArr_l[qttilde])-u*ndo_converged)) * (1.0/(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)-complex<double>(0.0,(2.0*wttilde)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkb]-kArr_l[qttilde])-u*(1.0-ndo_converged)));
                                                }
                                            }
                                            lower_level *= -1.0*SPINDEG*u/(beta*Nk); /// Removed minus sign
                                            //lower_level += 1.0;
                                            //cout << "gamma_oneD_spsp: " << u << " " << Nk << " " << lower_level << endl;
                                            return lower_level;//_Gk._u/lower_level
                                        };
                                        auto gamma_lambda_val = gamma_lambda_expr();
                                        complex<double> weights_val = (1.0/(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkt])-u*ndo_converged)
                                        )*(1.0/(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)+qq._iwn+u/2.0+2.0*cos(kArr_l[lkt]+qq._qx)-u*ndo_converged)
                                        )*(1.0/(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)+qq._iwn+u/2.0+2.0*cos(kArr_l[lkb]+qq._qx)-u*(1.0-ndo_converged))
                                        )*(1.0/(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)+u/2.0+2.0*cos(kArr_l[lkb])-u*(1.0-ndo_converged))
                                        );
                                        //lock_guard<mutex> guard(mutx);
                                        tmp_val_kt_kb += gamma_lambda_val;
                                        tmp_val_weigths += weights_val;
                                        if ((wtilde==0) && (wbar==0)){
                                            cout << "Thread id: " << this_thread::get_id() << endl;
                                            cout << "ktilde: " << kArr_l[lkt] << endl;
                                            cout << "kbar: " << kArr_l[lkb] << endl;
                                            cout << "gamma_oneD_spsp: " << tmp_val_kt_kb << endl;
                                            cout << "weights: " << tmp_val_weigths << endl;
                                        }
                                    } 
                                }
                                cout << "thread id operating: " << this_thread::get_id() << endl;
                            }); // END OF THREAD t1!!
                            tt[k]=move(t1);
                        }
                        join_all(tt);
                    }
                }
                it+=NUM_THREADS;
            }
            // arma::cx_dmat::iterator matGammaPtr = matGamma.begin();
            // arma::cx_dmat::iterator matWeigthsPtr = matWeigths.begin();
            // unsigned int it=0;
            // int totSize=kArr_l.size()*kArr_l.size();
            // cout << "totSize: " << totSize << " size kArr_l: " << kArr_l.size() << endl;
            // ThreadWrapper threadObj(u_ndo_c,qq);
            // // ThreadWrapper threadObj(u_ndo_c,qq,matGamma,matWeigths);
            // while (it<totSize){
            //     if (totSize % NUM_THREADS != 0){
            //         if ( totSize-it == (totSize % NUM_THREADS) ){
            //             int newl=totSize-it;
            //             vector<thread> tt(newl);
            //             for (int l=0; l<newl; l++){
            //                 int ltot=it+l; // Have to make sure spans over the whole array of k-space.
            //                 int lkt = static_cast<int>(floor(ltot/kArr_l.size()));
            //                 int lkb = (totSize % kArr_l.size());
            //                 thread t(ref(threadObj),lkt,lkb,beta);
            //                 tt[l]=move(t);
            //                 // tt[l]=thread(threadObj,lkt,lkb,beta);
            //                 cout << "ho" << ltot << "\n";
            //             }
            //             threadObj.join_all(tt);
            //         }
            //         else{
            //             vector<thread> tt(NUM_THREADS);
            //             for (int l=0; l<NUM_THREADS; l++){
            //                 int ltot=it+l; // Have to make sure spans over the whole array of k-space.
            //                 int lkt = static_cast<int>(floor(ltot/kArr_l.size()));
            //                 int lkb = (ltot % kArr_l.size());
            //                 cout << "lkt: " << lkt << " lkb: " << lkb << "\n";
            //                 thread t(ref(threadObj),lkt,lkb,beta);
            //                 tt[l]=move(t);
            //                 //tt.push_back(static_cast<thread&&>(t));
            //                 // tt[l]=thread(threadObj,lkt,lkb,beta);
            //                 cout << "hola" << ltot << "\n";
            //             }
            //             threadObj.join_all(tt);
            //         }
            //     }
            //     else{
            //         vector<thread> tt(NUM_THREADS);
            //         for (int l=0; l<NUM_THREADS; l++){
            //             int ltot=it+l; // Have to make sure spans over the whole array of k-space.
            //             int lkt = static_cast<int>(floor(ltot/kArr_l.size()));
            //             int lkb = (ltot % kArr_l.size());
            //             cout << "lkt: " << lkt << " lkb: " << lkb << "\n";
            //             thread t(ref(threadObj),lkt,lkb,beta);
            //             tt[l]=move(t);
            //             // tt[l]=thread(threadObj,lkt,lkb,beta);
            //         }
            //         threadObj.join_all(tt);
            //     }
            //     it+=NUM_THREADS;
            //     if (it > 6){ // To test speeds of different configurations...
            //         exit(0);
            //     }
            // }
            // outputFileChispspGamma.open(fileOutputChispspGamma, ofstream::out | ofstream::app);
            // outputFileChispspWeights.open(fileOutputChispspWeigths, ofstream::out | ofstream::app);
            // for (int ktilde=0; ktilde<kArr_l.size(); ktilde++){
            //     for (int kbar=0; kbar<kArr_l.size(); kbar++){
            //         outputFileChispspGamma << matGamma(kbar,ktilde) << " ";
            //         outputFileChispspWeights << matWeigths(kbar,ktilde) << " ";
            //     }
            //     outputFileChispspGamma << "\n";
            //     outputFileChispspWeights << "\n";
            // }
            // outputFileChispspGamma.close();
            // outputFileChispspWeights.close();
            #else
            Susceptibility susObj;
            for (int ktilde=0; ktilde<kArr_l.size(); ktilde++){
                cout << "ktilde: " << ktilde << endl;
                for (int kbar=0; kbar<kArr_l.size(); kbar++){
                    cout << "kbar: " << kbar << endl;
                    complex<double> tmp_val_kt_kb(0.0,0.0), tmp_val_kt_kb_bubble(0.0,0.0);
                    complex<double> tmp_val_weigths(0.0,0.0);
                    for (int wtilde=0; wtilde<Gup_k.size(); wtilde++){
                        for (int wbar=0; wbar<Gup_k.size(); wbar++){
                            tuple< complex<double>, complex<double> > gammaStuff;
                            if ( (wtilde==0) && (wbar==0) ){ // setting some conditions for the Matsubara frequencies (lowest frequencies and weights modified (beta)). same k grid plotted!
                                // cout << wtilde << " and wbar " << wbar << endl;
                                if (!is_full){
                                    gammaStuff=susObj.gamma_oneD_spsp_plotting(u_ndo_c,kArr_l[ktilde],complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta),kArr_l[kbar],complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta),qq);
                                    //gammaStuff=susObj.gamma_oneD_spsp_crossed_plotting(u_ndo_c,kArr_l[ktilde],complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta),kArr_l[kbar],complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta),qq);
                                }
                                else{
                                    gammaStuff=susObj.gamma_oneD_spsp_full_middle_plotting(u_ndo_c,kArr_l[kbar],kArr_l[ktilde],complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta),complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta),qq);
                                }
                                tmp_val_kt_kb += get<0>(gammaStuff);
                                tmp_val_kt_kb_bubble += get<1>(gammaStuff);
                                tmp_val_weigths += u_ndo_c(
                                    complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta),kArr_l[ktilde]
                                    )(0,0)*u_ndo_c(
                                    complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)+qq._iwn,kArr_l[ktilde]+qq._qx
                                    )(0,0)*u_ndo_c(
                                    complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)+qq._iwn,kArr_l[kbar]+qq._qx
                                    )(1,1)*u_ndo_c(
                                    complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta),kArr_l[kbar]
                                    )(1,1);
                            }
                        } 
                    }
                    // tmp_val_weigths *= 1.0/(beta)/(beta);
                    // tmp_val_kt_kb *= 1.0/(beta)/(beta);
                    // tmp_val_kt_kb_bubble *= 1.0/(beta)/(beta);
                    outputFileChispspGamma.open(fileOutputChispspGamma, ofstream::out | ofstream::app);
                    outputFileChispspWeights.open(fileOutputChispspWeigths, ofstream::out | ofstream::app);
                    outputFileChispspGammaBubble.open(fileOutputChispspGammaBubble, ofstream::out | ofstream::app);
                    outputFileChispspGamma << tmp_val_kt_kb << " ";
                    outputFileChispspGammaBubble << tmp_val_kt_kb_bubble << " ";
                    outputFileChispspWeights << tmp_val_weigths << " ";
                    outputFileChispspGamma.close();
                    outputFileChispspGammaBubble.close();
                    outputFileChispspWeights.close();
                }
                outputFileChispspGamma.open(fileOutputChispspGamma, ofstream::out | ofstream::app);
                outputFileChispspGammaBubble.open(fileOutputChispspGammaBubble, ofstream::out | ofstream::app);
                outputFileChispspWeights.open(fileOutputChispspWeigths, ofstream::out | ofstream::app);
                outputFileChispspGamma << "\n";
                outputFileChispspGammaBubble << "\n";
                outputFileChispspWeights << "\n";
                outputFileChispspGamma.close();
                outputFileChispspGammaBubble.close();
                outputFileChispspWeights.close();

            }
            
            #endif
            

            #else
            // complex<double> testSup;
            if (VERBOSE > 0) cout << "First 2D: " << u_ndo_c << endl;
            u_ndo_c.get_ndo_2D();
            cout << "After 2D: " << u_ndo_c << endl;

            Hubbard::K_2D qq(0.0,0.0,0.0+0.0*im); // photon 4-vector

            Susceptibility susObj;
            for (int kbary_m_tildey=0; kbary_m_tildey<kArr_l.size(); kbary_m_tildey++){
                cout << "ktildey_m_bary: " << kbary_m_tildey << "\n";
                for (int kbarx_m_tildex=0; kbarx_m_tildex<kArr_l.size(); kbarx_m_tildex++){
                    cout << "ktildex_m_barx: " << kbarx_m_tildex << "\n";
                    complex<double> tmp_val_kt_kb(0.0,0.0), tmp_val_kt_kb_bubble(0.0,0.0);
                    for (int wtilde=0; wtilde<Gup_k.size(); wtilde++){
                        for (int wbar=0; wbar<Gup_k.size(); wbar++){
                            tuple< complex<double>, complex<double> > gammaStuff;
                            if ( (wtilde==0) && (wbar==0) ){ // setting some conditions for the Matsubara frequencies (lowest frequencies and weights modified (beta)). same k grid plotted!
                                // cout << wtilde << " and wbar " << wbar << endl;
                                gammaStuff=susObj.gamma_twoD_spsp_plotting(u_ndo_c,kArr_l[kbarx_m_tildex],kArr_l[kbary_m_tildey],complex<double>(0.0,(2.0*(double)wtilde+1.0)*M_PI/beta),complex<double>(0.0,(2.0*(double)wbar+1.0)*M_PI/beta));
                                tmp_val_kt_kb += get<0>(gammaStuff);
                                tmp_val_kt_kb_bubble += get<1>(gammaStuff);
                            }
                        }
                    }
                    // tmp_val_weigths *= 1.0/(beta)/(beta);
                    // tmp_val_kt_kb *= 1.0/(beta)/(beta);
                    // tmp_val_kt_kb_bubble *= 1.0/(beta)/(beta);
                    outputFileChispspGamma.open(fileOutputChispspGamma, ofstream::out | ofstream::app);
                    outputFileChispspGammaBubble.open(fileOutputChispspGammaBubble, ofstream::out | ofstream::app);
                    outputFileChispspGamma << tmp_val_kt_kb << " ";
                    outputFileChispspGammaBubble << tmp_val_kt_kb_bubble << " ";
                    outputFileChispspGamma.close();
                    outputFileChispspGammaBubble.close();
                }
                outputFileChispspGamma.open(fileOutputChispspGamma, ofstream::out | ofstream::app);
                outputFileChispspGammaBubble.open(fileOutputChispspGammaBubble, ofstream::out | ofstream::app);
                outputFileChispspGamma << "\n";
                outputFileChispspGammaBubble << "\n";
                outputFileChispspGamma.close();
                outputFileChispspGammaBubble.close();
            }
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


inline bool file_exists (const string& filename) {
  struct stat buffer;   
  return (stat(filename.c_str(), &buffer) == 0); 
}

void getWeights(Hubbard::FunctorBuildGk& Gk, Hubbard::K_2D qq, vector<double>& kArr_l, double beta, double Nomega, ofstream&outputFileChispspWeights, string& fileOutputChispspWeigths){
    for (int ktildey=0; ktildey<kArr_l.size(); ktildey++){
        for (int kbary=0; kbary<kArr_l.size(); kbary++){
            for (int ktildex=0; ktildex<kArr_l.size(); ktildex++){
                for (int kbarx=0; kbarx<kArr_l.size(); kbarx++){
                    complex<double> tmp_val_weigths(0.0,0.0);
                    for (int wtilde=0; wtilde<Nomega; wtilde++){
                        for (int wbar=0; wbar<Nomega; wbar++){
                            if ( (wtilde==0) || (wbar==0) ){
                                tmp_val_weigths += Gk(
                                    complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta),kArr_l[ktildex],kArr_l[ktildey]
                                    )(0,0)*Gk(
                                    complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)+qq._iwn,kArr_l[ktildex]+qq._qx,kArr_l[ktildey]+qq._qy
                                    )(0,0)*Gk(
                                    complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)+qq._iwn,kArr_l[kbarx]+qq._qx,kArr_l[kbary]+qq._qy
                                    )(1,1)*Gk(
                                    complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta),kArr_l[kbarx],kArr_l[kbary]
                                    )(1,1);
                            }
                        }
                    }
                    outputFileChispspWeights.open(fileOutputChispspWeigths, ofstream::out | ofstream::app);
                    outputFileChispspWeights << tmp_val_weigths << " ";
                    outputFileChispspWeights.close();
                }
            }
        }
    }
    outputFileChispspWeights.open(fileOutputChispspWeigths, ofstream::out | ofstream::app);
    outputFileChispspWeights << "\n";
    outputFileChispspWeights.close();
}

#ifdef PARALLEL

void join_all(vector<thread>& grp){
    for (auto& thread : grp){
        if (thread.joinable()){
            thread.join();
        }
    }
}

// ThreadWrapper::ThreadWrapper(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D& q,arma::cx_dmat::iterator matPtr,arma::cx_dmat::iterator matWPtr) : 
// _Gk(Gk), _q(q), _ktb(matPtr), _ktbW(matWPtr) {}

// ThreadWrapper::ThreadWrapper(Hubbard::FunctorBuildGk& Gk, Hubbard::K_1D& q) : 
// _Gk(Gk), _q(q), _ktb(nullptr), _ktbW(nullptr) {}

// void ThreadWrapper::operator()(int ktilde, int kbar, double beta){
//     complex<double> tmp_val_kt_kb(0.0,0.0);
//     complex<double> tmp_val_weigths(0.0,0.0);
//     // cout << beta << " " << _Gk._kArr_l[ktilde] << " " << _q._iwn << " " << _q._qx << " " << _Gk._kArr_l[kbar] << endl;
//     for (int wtilde=0; wtilde<_Gk._size; wtilde++){
//         for (int wbar=0; wbar<_Gk._size; wbar++){
//             //lock_guard<mutex> guard(mutx);
//             //tmp_val_kt_kb += gamma_oneD_spsp(_Gk._kArr_l[ktilde],complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta),_Gk._kArr_l[kbar],complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta));
//             auto gamma_lambda_expr = [&](){
//                 complex<double> lower_level=0.0+0.0*im;
//                 for (int wttilde=0; wttilde<_Gk._size; wttilde++){
//                     for (size_t qttilde=0; qttilde<_Gk._kArr_l.size(); qttilde++){
//                         lower_level += _Gk(complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)-_Gk._precomp_qn[wttilde],_Gk._kArr_l[ktilde]-_Gk._kArr_l[qttilde])(0,0)*_Gk(complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)-_Gk._precomp_qn[wttilde],_Gk._kArr_l[kbar]-_Gk._kArr_l[qttilde])(1,1);
//                     }
//                 }
//                 lower_level *= -1.0*SPINDEG*_Gk._u/(_Gk._beta*_Gk._Nk); /// Removed minus sign
//                 //lower_level += 1.0;
//                 //cout << "gamma_oneD_spsp: " << _Gk._u << " " << _Gk._Nk << endl;
//                 return lower_level;//_Gk._u/lower_level
//             };
//             //lock_guard<mutex> guard(mutx);
//             tmp_val_kt_kb += gamma_lambda_expr();
//             tmp_val_weigths += _Gk(
//                                     complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta),_Gk._kArr_l[ktilde]
//                                     )(0,0)*_Gk(
//                                     complex<double>(0.0,(2.0*wtilde+1.0)*M_PI/beta)-_q._iwn,_Gk._kArr_l[ktilde]-_q._qx
//                                     )(0,0)*_Gk(
//                                     complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta)+_q._iwn,_Gk._kArr_l[kbar]+_q._qx
//                                     )(1,1)*_Gk(
//                                     complex<double>(0.0,(2.0*wbar+1.0)*M_PI/beta),_Gk._kArr_l[kbar]
//                                     )(1,1);
//             if ((wtilde==0) && (wbar==0)){
//                 cout << "Thread id: " << this_thread::get_id() << endl;
//                 cout << "ktilde: " << _Gk._kArr_l[ktilde] << endl;
//                 cout << "kbar: " << _Gk._kArr_l[kbar] << endl;
//                 cout << "gamma_oneD_spsp: " << tmp_val_kt_kb << endl;
//                 cout << "weights: " << tmp_val_weigths << endl;
//             }
//         } 
//     }
//     // lock_guard<mutex> guard(mutx);
//     matGamma(kbar,ktilde) = tmp_val_kt_kb; // These matrices are static variables.
//     matWeigths(kbar,ktilde) = 1.0/(_Gk._beta*_Gk._Nk)/(_Gk._beta*_Gk._Nk)*tmp_val_weigths;
//     //cout << "Gamma for " << "ktilde " << ktilde << " and kbar " << kbar << ": " << matGamma(kbar,ktilde) << "\n";
//     //cout << "Weigths for " << "ktilde " << ktilde << " and kbar " << kbar << ": " << matWeigths(kbar,ktilde) << "\n";
// }

// complex<double> ThreadWrapper::gamma_oneD_spsp(double ktilde,complex<double> wtilde,double kbar,complex<double> wbar){
//     complex<double> lower_level=0.0+0.0*im;
//     for (int wttilde=0; wttilde<_Gk._size; wttilde++){
//         for (size_t qttilde=0; qttilde<_Gk._kArr_l.size(); qttilde++){
//             lower_level += _Gk(wtilde-_Gk._precomp_qn[wttilde],ktilde-_Gk._kArr_l[qttilde])(0,0)*_Gk(wbar-_Gk._precomp_qn[wttilde],kbar-_Gk._kArr_l[qttilde])(1,1);
//         }
//     }
    
//     lower_level *= -1.0*SPINDEG*_Gk._u/(_Gk._beta*_Gk._Nk); /// Removed minus sign
//     //lower_level += 1.0;
//     //cout << "gamma_oneD_spsp: " << _Gk._u/lower_level << endl;
//     return lower_level;//_Gk._u/lower_level;
// }


// void ThreadWrapper::join_all(vector<thread>& grp){
//     for (auto& thread : grp){
//         if (thread.joinable()){
//             thread.join();
//         }
//     }
// }

#endif
