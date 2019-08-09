#include "susceptibility_utils.h"


std::complex<double> Susceptibility::gamma_oneD_spsp(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D ktilde,Hubbard::K_1D kbar,Hubbard::K_1D q) const{ // q's contain bosonic Matsubara frequencies.
    std::complex<double> lower_level=0.0+0.0*im;
    for (int wttilde=0; wttilde<Gk._size; wttilde++){
        for (size_t qttilde=0; qttilde<Gk._kArr_l.size(); qttilde++){
            lower_level += Gk(ktilde._iwn-Gk._precomp_qn[wttilde],ktilde._qx-Gk._kArr_l[qttilde])(0,0)*Gk((kbar+q)._iwn-Gk._precomp_qn[wttilde],(kbar+q)._qx-Gk._kArr_l[qttilde])(1,1);
        }
    }
    lower_level *= -1.0*Gk._u/(Gk._beta*Gk._Nk); /// Removed minus sign
    lower_level += 1.0;
    return Gk._u/lower_level;
}

std::complex<double> Susceptibility::chispsp(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const{
    std::complex<double> upper_level=0.0+0.0*im;
    std::ofstream output("ktilde_kbar_U"+std::to_string(Gk._u)+"_beta"+std::to_string(Gk._beta)+"_ndo"+std::to_string(Gk._ndo)+"_Nk"+std::to_string(Gk._Nk)+"_Nw"+std::to_string(Gk._size)+".dat");
    for (size_t ktilde=0; ktilde<Gk._kArr_l.size(); ktilde++){
        std::cout << "ktilde: " << ktilde << std::endl;
        for (size_t kbar=0; kbar<Gk._kArr_l.size(); kbar++){
            std::cout << "kbar: " << kbar << std::endl;
            std::complex<double> tmp_val_kbar=0.0+0.0*im;
            for (int wtilde=0; wtilde<Gk._size; wtilde++){
                for (int wbar=0; wbar<Gk._size; wbar++){
                    Hubbard::K_1D kt(Gk._kArr_l[ktilde], Gk._precomp_wn[wtilde]);
                    Hubbard::K_1D kb(Gk._kArr_l[kbar],Gk._precomp_wn[wbar]);
                    std::complex<double> tmp_val = Gk(
                            Gk._precomp_wn[wtilde],Gk._kArr_l[ktilde]
                            )(0,0)*Gk(
                            Gk._precomp_wn[wtilde]-q._iwn,Gk._kArr_l[ktilde]-q._qx
                            )(0,0)*gamma_oneD_spsp( Gk, kt, kb, q )*Gk(
                            Gk._precomp_wn[wbar]+q._iwn,Gk._kArr_l[kbar]+q._qx
                            )(1,1)*Gk(
                            Gk._precomp_wn[wbar],Gk._kArr_l[kbar]
                            )(1,1);
                    upper_level += tmp_val;
                    tmp_val_kbar += tmp_val;
                }
            }
            output << abs(tmp_val_kbar) << " ";
        }
        output << "\n";
    }
    output.close();
    upper_level *= -1.0*(1.0/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk));
    return upper_level;
}

std::complex<double> Susceptibility::gamma_oneD_spsp(Hubbard::FunctorBuildGk& Gk,double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar,Hubbard::K_1D q) const{ // q's contain bosonic Matsubara frequencies.
    std::complex<double> lower_level=0.0+0.0*im;
    for (int wttilde=0; wttilde<Gk._size; wttilde++){
        for (size_t qttilde=0; qttilde<Gk._kArr_l.size(); qttilde++){
            lower_level += Gk(wtilde-Gk._precomp_qn[wttilde],ktilde-Gk._kArr_l[qttilde])(0,0)*Gk(wbar+q._iwn-Gk._precomp_qn[wttilde],kbar+q._qx-Gk._kArr_l[qttilde])(1,1);
        }
    }
    lower_level *= -1.0*Gk._u/(Gk._beta*Gk._Nk); /// Removed minus sign
    lower_level += 1.0;
    return Gk._u/lower_level;
}

std::complex<double> Susceptibility::chispsp_long_expr(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const{
    std::ofstream output;
    std::string strOutput("ktilde_kbar_U"+std::to_string(Gk._u)+"_beta"+std::to_string(Gk._beta)+"_ndo"+std::to_string(Gk._ndo)+"_Nk"+std::to_string(Gk._Nk)+"_Nw"+std::to_string(Gk._size)+"__AA.dat");
    std::complex<double> upper_level=0.0+0.0*im;
    for (size_t ktilde=0; ktilde<Gk._kArr_l.size(); ktilde++){
        std::cout << "ktilde: " << ktilde << std::endl;
        for (size_t kbar=0; kbar<Gk._kArr_l.size(); kbar++){
            std::cout << "kbar: " << kbar << std::endl;
            std::complex<double> tmp_val_kbar=0.0+0.0*im;
            for (int wtilde=0; wtilde<Gk._size; wtilde++){
                for (int wbar=0; wbar<Gk._size; wbar++){
                    std::complex<double> tmp_val = Gk(
                            Gk._precomp_wn[wtilde],Gk._kArr_l[ktilde]
                            )(0,0)*Gk(
                            Gk._precomp_wn[wtilde]-q._iwn,Gk._kArr_l[ktilde]-q._qx
                            )(0,0)*gamma_oneD_spsp( Gk, Gk._kArr_l[ktilde], Gk._precomp_wn[wtilde], Gk._kArr_l[kbar], Gk._precomp_wn[wbar], q )*Gk(
                            Gk._precomp_wn[wbar]+q._iwn,Gk._kArr_l[kbar]+q._qx
                            )(1,1)*Gk(
                            Gk._precomp_wn[wbar],Gk._kArr_l[kbar]
                            )(1,1);
                    upper_level += tmp_val;
                    tmp_val_kbar += tmp_val;
                }
                
            }
            output.open(strOutput, std::ofstream::out | std::ofstream::app);
            output << tmp_val_kbar << " ";
            output.close();
        }
        output.open(strOutput, std::ofstream::out | std::ofstream::app);
        output << "\n";
        output.close();
    }
    upper_level *= -1.0*(1.0/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk));
    output.open(strOutput, std::ofstream::out | std::ofstream::app);
    output << "\n\n";
    output << "Spin susceptibility at "+std::to_string(Gk._u)+" : " << upper_level << "\n\n";
    output.close();
    return upper_level;
}

std::complex<double> Susceptibility::chi0(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const{
    std::complex<double> xi=0.0+0.0*im;
    for (size_t ikn=0; ikn<Gk._size; ikn++){
        for (size_t k=0; k<Gk._kArr_l.size(); k++){
            xi += Gk( Gk._precomp_wn[ikn] + q._iwn, Gk._kArr_l[k] + q._qx )(0,0) * Gk( Gk._precomp_wn[ikn], Gk._kArr_l[k] )(0,0);
        }
    }
    xi *= 1.0/(Gk._beta*Gk._Nk); /// Removed minus sign
    return xi;
}
        
std::complex<double> Susceptibility::chi0(Hubbard::FunctorBuildGk& Gk,Hubbard::K_2D q) const{
    std::complex<double> xi=0.0+0.0*im;
    for (size_t ikn=0; ikn<Gk._size; ikn++){
        for (size_t kx=0; kx<Gk._kArr_l.size(); kx++){
            for (size_t ky=0; ky<Gk._kArr_l.size(); ky++){
                xi += Gk(Gk._precomp_wn[ikn]+q._iwn,Gk._kArr_l[kx]+q._qx,Gk._kArr_l[ky]+q._qy)(0,0)*Gk(Gk._precomp_wn[ikn],Gk._kArr_l[kx],Gk._kArr_l[ky])(0,0);
            }
        }
    }
    xi *= 1.0/(Gk._beta*Gk._Nk*Gk._Nk); /// Removed minus sign
    return xi;
}

std::complex<double> Susceptibility::chisp(Hubbard::FunctorBuildGk& Gk, Hubbard::K_1D q) const{
    std::complex<double> xi=0.0+0.0*im, suscept(0.0,0.0);
    xi=chi0(Gk,q);
    suscept = xi/(1.0-Gk._u*xi);
    return suscept;
}

std::complex<double> Susceptibility::chisp(Hubbard::FunctorBuildGk& Gk, Hubbard::K_2D q) const{
    std::complex<double> xi=0.0+0.0*im, suscept(0.0,0.0);
    xi=chi0(Gk,q); /// Removed minus sign
    suscept = xi/(1.0-Gk._u*xi);
    return suscept;
}