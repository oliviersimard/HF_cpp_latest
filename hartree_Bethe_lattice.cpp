#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <string>

void init_G0(std::vector<std::complex<double> >&,int,int,double,double);
void solv_dyson_eq_up(std::vector< std::complex<double> >&,std::vector< std::complex<double> >&,double,double);
void solv_dyson_eq_do(std::vector< std::complex<double> >&,std::vector< std::complex<double> >&,double,double);
void latt_dyson_eq_up(std::vector< std::complex<double> >&,std::vector< std::complex<double> >&,double,double);
void latt_dyson_eq_do(std::vector< std::complex<double> >&,std::vector< std::complex<double> >&,double,double);

int main(int argc, char ** argv){
 
	double u,beta,ndo_av;
	int Nomega=400;
    double ndo_initial=0.6;
    int Niterations=100;
    double mu=0.0;
    int N_e=1000; // Energy discretization
    double ndo=ndo_initial;
    std::vector<std::complex<double> > Gup_k(Nomega), Gdo_k(Nomega);
    std::vector<std::complex<double> > Gup_k_weiss(Nomega), Gdo_k_weiss(Nomega);
  
    std::ofstream output;
    std::string strOutput("beta_5.0_200.0_1.0_vs_u_0.0_3.0_0.01_Bethe_lattice.dat");
    for (beta=5.0; beta<=200.0; beta+=1.0){
        
        for (u=0.0; u<=3.0; u+=0.01) {
            mu=u/2.;
            ndo=ndo_initial;
            // Initialize the weiss Green's function with semicircular DOS. (Bethe lattice)
            init_G0(Gup_k_weiss,Nomega,N_e,beta,mu);
            init_G0(Gdo_k_weiss,Nomega,N_e,beta,mu);
    
            for (int i=0; i<Niterations; i++) {
            
                ndo_av=0;

                // Get dressed Gup_k
                solv_dyson_eq_up(Gup_k_weiss,Gup_k,ndo,u);
                // Get dressed Gdo_k
                solv_dyson_eq_do(Gdo_k_weiss,Gdo_k,ndo,u);

                for (int j=0; j<Gup_k.size(); j++)
                    ndo_av += (2./beta)*(Gup_k[j]-1./std::complex<double>(0,(2*j+1)*M_PI/beta)).real();
                ndo_av -= 0.5;
                ndo_av *= (-1);

                // Update the Weiss Green's functions
                latt_dyson_eq_up(Gup_k_weiss, Gdo_k, mu, beta);
                latt_dyson_eq_do(Gdo_k_weiss, Gup_k, mu, beta);
            
                ndo = ndo_av; // Updated ndo
            }
            output.open(strOutput, std::ofstream::out | std::ofstream::app);
            output << ndo << " ";
            output.close();
            std::cout << "beta: " << beta << " u: " << u << " ndo: " << ndo << "\n";
        }
        std::cout << "\n";
        output.open(strOutput, std::ofstream::out | std::ofstream::app);
        output << "\n";
        output.close();
    }
  return 0;
}


void init_G0(std::vector<std::complex<double> >& Gup_k_weiss, int Nomega, int N_e, double beta, double mu){
    double e_max=2.0;
    double de=e_max/(double)N_e;
    for (int j=0; j<Nomega; j++){
        Gup_k_weiss[j] = std::complex<double>(0.0,0.0);
        std::complex<double> w = std::complex<double>(0.0,(2.0*(double)j+1.0)*M_PI/beta);
        for (int k=1; k<2*N_e; k++){
            double ek=(k-N_e)*de;
            Gup_k_weiss[j] += de*sqrt(4.0-ek*ek)/(2.0*M_PI)/(w+mu-ek);
        }
    }
}

void solv_dyson_eq_up(std::vector< std::complex<double> >& Gup_k_weiss, std::vector< std::complex<double> >& Gup_k, double ndo, double U){
    for (size_t j=0; j<Gup_k.size(); j++){
        Gup_k[j] = Gup_k_weiss[j]/(1.0-Gup_k_weiss[j]*U*ndo);
    }
}

void solv_dyson_eq_do(std::vector< std::complex<double> >& Gdo_k_weiss, std::vector< std::complex<double> >& Gdo_k, double ndo, double U){
    for (size_t j=0; j<Gdo_k.size(); j++){
        Gdo_k[j] = Gdo_k_weiss[j]/(1.0-Gdo_k_weiss[j]*U*(1.0-ndo));
    }
}

void latt_dyson_eq_up(std::vector< std::complex<double> >& Gup_k_weiss, std::vector< std::complex<double> >& Gdo_k, double mu, double beta){
    for (size_t j=0; j<Gup_k_weiss.size(); j++){
        Gup_k_weiss[j] = 1.0/( std::complex<double>(0.0,(2.0*(double)j+1.0)*M_PI/beta) + mu - Gdo_k[j] );
    }
}

void latt_dyson_eq_do(std::vector< std::complex<double> >& Gdo_k_weiss, std::vector< std::complex<double> >& Gup_k, double mu, double beta){
    for (size_t j=0; j<Gdo_k_weiss.size(); j++){
        Gdo_k_weiss[j] = 1.0/( std::complex<double>(0.0,(2.0*(double)j+1.0)*M_PI/beta) + mu - Gup_k[j] );
    }
}
