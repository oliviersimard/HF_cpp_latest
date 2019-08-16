#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <string>
#include <armadillo>

const int _Nomega=100,_Niterations=200,_Nk=100;

typedef arma::Mat< std::complex<double> > arma_t;

const arma_t ZEROS = arma::Mat< std::complex<double> >(_Nomega,_Nk,arma::fill::zeros);

class HubbardExt{
    public:
        HubbardExt(double,double,double,double,int,int,int);
        //void operator()(int,std::complex<double>);
        void init_self();
        void Gup_kAA(int,int,int,int);
        void Sup_kAA(int,int);
        void Gdo_kAA(int,int,int,int);
        void Sdo_kAA(int,int);
        void Gup_kBB(int,int,int,int);
        void Sup_kBB(int,int);
        void Gdo_kBB(int,int,int,int);
        void Sdo_kBB(int,int);
        std::complex<double> w(int);
        double get_nAA();
        double get_double_occupancy_AA();
    private:
        double _u,_v,_beta,_mu;
        std::vector<double> _kArr;
        arma_t _Gup_kAA, _Gdo_kAA, _Gup_kBB, _Gdo_kBB;
        arma_t _Sup_kAA, _Sdo_kAA, _Sup_kBB, _Sdo_kBB;
        std::vector< std::complex<double>* > _G_vec_ptr;
        std::vector< std::complex<double>* > _S_vec_ptr;
        // Static members
        static std::complex<double> self_init;
};

int main(int argc, char ** argv){
 
	double n,d;
    double beta_init=5.0, beta_step=5.0, beta_max=100.0;
    double u_init=2.0, u_step=0.2, u_max=4.0;
    double v=0.5;
    double mu=0.0;
  
    std::ofstream output;
    std::ofstream outputDO;
    std::string strOutput("T_vs_U_beta_"+std::to_string(beta_init)+"_"+std::to_string(beta_max)+"_"+std::to_string(beta_step)+"_vs_u_"+std::to_string(u_init)+"_"+std::to_string(u_max)+"_"+std::to_string(u_step)+"_v_"+std::to_string(v)+"_Nk_"+std::to_string(_Nk)+"_Nomega_"+std::to_string(_Nomega)+"_Nit_"+std::to_string(_Niterations)+"_Fock.dat");
    std::string strOutputDO("DO_beta_"+std::to_string(beta_init)+"_"+std::to_string(beta_max)+"_"+std::to_string(beta_step)+"_vs_u_"+std::to_string(u_init)+"_"+std::to_string(u_max)+"_"+std::to_string(u_step)+"_v_"+std::to_string(v)+"_Nk_"+std::to_string(_Nk)+"_Nomega_"+std::to_string(_Nomega)+"_Nit_"+std::to_string(_Niterations)+"_Fock.dat");
    for (double beta=beta_init; beta<=beta_max; beta+=beta_step){
        
        for (double u=u_init; u<=u_max; u+=u_step) {
            mu=u/2.0+2.0*v; // chemical potential in the 1D case for the extended Hubbard model.

            HubbardExt hubbardExtObj(u,v,beta,mu,_Nomega,_Niterations,_Nk);
            // Initialize the self-energies.
            hubbardExtObj.init_self();
    
            for (int i=0; i<_Niterations; i++){
                
                for (int kn=0; kn<_Nomega; kn++){
                    
                    for (int k=0; k<_Nk; k++){
                        // First updating the AA self-energies
                        hubbardExtObj.Sup_kAA(k,kn);
                        hubbardExtObj.Sdo_kAA(k,kn);
                    }
                }
                for (int kn=0; kn<_Nomega; kn++){
                    
                    for (int k=0; k<_Nk; k++){
                        // Then updating the BB self-energies
                        hubbardExtObj.Sup_kBB(k,kn);
                        hubbardExtObj.Sdo_kBB(k,kn);
                    }
                }
                std::cout << "it: " << i << std::endl;
                n = hubbardExtObj.get_nAA();
                d = hubbardExtObj.get_double_occupancy_AA();
                std::cout << "d: " << d << std::endl;
            }
            output.open(strOutput, std::ofstream::out | std::ofstream::app);
            output << n << " ";
            output.close();

            outputDO.open(strOutputDO, std::ofstream::out | std::ofstream::app);
            outputDO << d << " ";
            outputDO.close();
        }
        std::cout << "\n";
        output.open(strOutput, std::ofstream::out | std::ofstream::app);
        output << "\n";
        output.close();

        outputDO.open(strOutputDO, std::ofstream::out | std::ofstream::app);
        outputDO << " ";
        outputDO.close();
    }
  return 0;
}

std::complex<double> HubbardExt::self_init = std::complex<double>(0.2,-0.1);

HubbardExt::HubbardExt(double u, double v, double beta, double mu, int Nomega, int Nit, int Nk) : _u(u), _v(v), _beta(beta), _mu(mu){
    for (int k=0; k<=_Nk; k++){
        this->_kArr.push_back(-M_PI/2.0 + 2.0*(double)k*M_PI/(2.0*(double)Nk));
    }
    _Gup_kAA=ZEROS; _Gdo_kAA=ZEROS; _Gup_kBB=ZEROS; _Gdo_kBB=ZEROS;
    _Sup_kAA=ZEROS; _Sdo_kAA=ZEROS; _Sup_kBB=ZEROS; _Sdo_kBB=ZEROS;
    this->_G_vec_ptr = { _Gup_kAA.memptr(), _Gdo_kAA.memptr(), _Gup_kBB.memptr(), _Gdo_kBB.memptr() };
    this->_S_vec_ptr = { _Sup_kAA.memptr(), _Sdo_kAA.memptr(), _Sup_kBB.memptr(), _Sdo_kBB.memptr() };
}

void HubbardExt::init_self(){
    for (size_t it=0; it<_S_vec_ptr.size(); it++){
        for (int j=0; j<_Nomega; j++){
            std::complex<double> w = std::complex<double>(0.0,(2.0*(double)j+1.0)*M_PI/_beta);
            for (int k=0; k<_kArr.size(); k++){
                *(_S_vec_ptr[it] + j*_Nk + k) = std::complex<double>(0.4,-0.1);// /(w+_mu-self_init);
            }
        }
        self_init += std::complex<double>(0.1,0.1);
    }
}

void HubbardExt::Gup_kAA(int k, int kn, int q, int qn){
    double epsk=-2.0*cos(_kArr[k]+_kArr[q]);
    std::complex<double> wj = w(kn)+w(kn);
    // if (Nit != 0)
    //     Sup_kAA(k,kn,); // If Nit != 
    *(_G_vec_ptr[0] + qn*_Nk + q) = 1.0/( wj + _mu - *(_S_vec_ptr[0] + qn*_Nk + q) - epsk*epsk/( wj + _mu - *(_S_vec_ptr[2] + qn*_Nk + q) ) );
}

void HubbardExt::Sup_kAA(int k, int kn){
    for (int qn=0; qn<_Nomega; qn++){
        for (int q=0; q<_Nk; q++){
            Gup_kAA(k,kn,q,qn);
            *(_S_vec_ptr[0] + kn*_Nk + k) += _v*cos(_kArr[q])*( -1.0/(_Nk*_beta) )*( *(_G_vec_ptr[0] + qn*_Nk + q) );
        }
    }
}

void HubbardExt::Gdo_kAA(int k, int kn, int q, int qn){
    double epsk=-2.0*cos(_kArr[k]+_kArr[q]);
    std::complex<double> wj = w(kn)+w(kn);
    *(_G_vec_ptr[1] + qn*_Nk + q) = 1.0/( wj + _mu - *(_S_vec_ptr[1] + qn*_Nk + q) - epsk*epsk/( wj + _mu - *(_S_vec_ptr[3] + qn*_Nk + q) ) );
}

void HubbardExt::Sdo_kAA(int k, int kn){
    for (int qn=0; qn<_Nomega; qn++){
        for (int q=0; q<_Nk; q++){
            Gdo_kAA(k,kn,q,qn);
            *(_S_vec_ptr[1] + kn*_Nk + kn) += _v*cos(_kArr[q])*( -1.0/(_Nk*_beta) )*( *(_G_vec_ptr[1] + qn*_Nk + q) );
        }
    }
}

// // BB quantities are computed after AA quantities have been computed.

void HubbardExt::Gup_kBB(int k, int kn, int q, int qn){
    double epsk=-2.0*cos(_kArr[k]+_kArr[q]);
    std::complex<double> wj = w(kn)+w(kn);
    *(_G_vec_ptr[2] + qn*_Nk + q) = 1.0/( wj + _mu - *(_S_vec_ptr[2] + qn*_Nk + q) - epsk*epsk/( wj + _mu - *(_S_vec_ptr[0] + qn*_Nk + q) ) );
}

void HubbardExt::Sup_kBB(int k, int kn){
    for (int qn=0; qn<_Nomega; qn++){
        for (int q=0; q<_Nk; q++){
            Gup_kBB(k,kn,q,qn);
            *(_S_vec_ptr[2] + kn*_Nk + k) += _v*cos(_kArr[q])*( -1.0/(_Nk*_beta) )*( *(_G_vec_ptr[2] + qn*_Nk + q) );
        }
    }
}

void HubbardExt::Gdo_kBB(int k, int kn, int q, int qn){
    double epsk=-2.0*cos(_kArr[k]+_kArr[q]);
    std::complex<double> wj = w(kn)+w(kn);
    *(_G_vec_ptr[3] + qn*_Nk + q) = 1.0/( wj + _mu - *(_S_vec_ptr[3] + qn*_Nk + q) - epsk*epsk/( wj + _mu - *(_S_vec_ptr[1] + qn*_Nk + q) ) );
}

void HubbardExt::Sdo_kBB(int k, int kn){
    for (int qn=0; qn<_Nomega; qn++){
        for (int q=0; q<_Nk; q++){
            Gup_kBB(k,kn,q,qn);
            *(_S_vec_ptr[3] + kn*_Nk + k) += _v*cos(_kArr[q])*( -1.0/(_Nk*_beta) )*( *(_G_vec_ptr[3] + qn*_Nk + q) );
        }
    }
}

// Once all the self-energies of one iteration have been calculated, build the new Green's functions out of the latter.

double HubbardExt::get_nAA(){
    double n=0.0;
    for (int j=0; j<_Nomega; j++){
        std::complex<double> n_k(0.0,0.0);
        std::complex<double> wj = w(j);
        for (int k=0; k<_kArr.size(); k++){ // Summing over G(k)
            double epsk = -2.0*cos(_kArr[k]);
            if ( (k==0) || (k==_Nk) ){
                n_k += 0.5/( wj + _mu - *(_S_vec_ptr[0] + j*_Nk + k) - epsk*epsk/( wj + _mu - *(_S_vec_ptr[2] + j*_Nk + k) ) );
            }
            else{
                n_k += 1.0/( wj + _mu - *(_S_vec_ptr[0] + j*_Nk + k) - epsk*epsk/( wj + _mu - *(_S_vec_ptr[2] + j*_Nk + k) ) );
            }
        }
        n_k /= (_Nk);
        n += (2.0/_beta)*( n_k - 1.0/wj ).real();
    }
    n -= 0.5;
    n *= -1.0;

    std::cout << "beta: " << _beta << " u: " << _u << " v: " << _v << " n: " << n << std::endl;

    return n;
}


double HubbardExt::get_double_occupancy_AA(){
/* This function computes U<n_{up}n_{down}> = 1/(\beta*V)*\sum_{k,ikn} \Sigma(k,ikn)G(k,ikn)e^{-ikn0^-} */
    double d=0.0;
    for (int j=0; j<_Nomega; j++){
        std::complex<double> d_k(0.0,0.0);
        std::complex<double> wj = w(j);
        for (int k=0; k<_kArr.size(); k++){ // Summing over G(k)
            double epsk = -2.0*cos(_kArr[k]);
            if ( (k==0) || (k==_Nk) ){
                d_k += *(_S_vec_ptr[0] + j*_Nk + k)*0.5/( wj + _mu - *(_S_vec_ptr[0] + j*_Nk + k) - epsk)//*epsk/( wj + _mu - *(_S_vec_ptr[2] + j*_Nk + k) ) );
            }
            else{
                d_k += *(_S_vec_ptr[0] + j*_Nk + k)*1.0/( wj + _mu - *(_S_vec_ptr[0] + j*_Nk + k) - epsk)//*epsk/( wj + _mu - *(_S_vec_ptr[2] + j*_Nk + k) ) );
            }
        }
        d_k /= (_Nk);
        d += (2.0/_beta)*( d_k ).real();
    }
    d *= -1.0;

    return d;
}

std::complex<double> HubbardExt::w(int j){
    return std::complex<double>(0.0,(2.0*(double)j+1.0)*M_PI/_beta);
}