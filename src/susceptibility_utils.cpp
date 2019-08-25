#include "susceptibility_utils.h"

/* Functions entering the full spin susceptibility. */

std::complex<double> Susceptibility::gamma_oneD_spsp_full_lower(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D ktilde,Hubbard::K_1D kbar,Hubbard::K_1D qtilde,Hubbard::K_1D q) const{
    std::complex<double> lower_level=0.0+0.0*im;
    for (int iqnttilde=0; iqnttilde<Gk._size; iqnttilde++){
        for (size_t qttilde=0; qttilde<Gk._kArr_l.size(); qttilde++){
            lower_level += Gk((ktilde+qtilde-q)._iwn,(ktilde+qtilde-q)._qx)(0,0)*Gk(kbar._iwn-Gk._precomp_qn[iqnttilde],kbar._qx-Gk._kArr_l[qttilde])(1,1);
        }
    }
    lower_level *= -1.0*Gk._u/(Gk._beta*Gk._Nk); /// Removed minus sign
    lower_level += 1.0;
    return Gk._u/lower_level;
}

std::complex<double> Susceptibility::gamma_oneD_spsp_full_middle(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D kbar,Hubbard::K_1D q) const{
    std::complex<double> middle_level=0.0+0.0*im;
    for (int iqntilde=0; iqntilde<Gk._size; iqntilde++){
        for (size_t qqtilde=0; qqtilde<Gk._kArr_l.size(); qqtilde++){
            for (int ikntilde=0; ikntilde<Gk._size; ikntilde++){
                for (size_t kktilde=0; kktilde<Gk._kArr_l.size(); kktilde++){
                    Hubbard::K_1D qtilde(Gk._kArr_l[qqtilde],Gk._precomp_qn[iqntilde]);
                    Hubbard::K_1D ktilde(Gk._kArr_l[kktilde],Gk._precomp_qn[ikntilde]);
                    middle_level += Gk(ktilde._iwn,ktilde._qx
                    )(0,0) * gamma_oneD_spsp_full_lower(Gk,ktilde,kbar,qtilde,q
                    ) * Gk(ktilde._iwn-q._iwn,ktilde._qx-ktilde._qx
                    )(0,0);
                }
            }
        }
    }
    middle_level *= -1.0*Gk._u/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk); /// Removed minus sign
    middle_level += 1.0;
    return middle_level;
}

std::complex<double> Susceptibility::chisp_full(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const{
    std::complex<double> upper_level=0.0+0.0*im;
    for (int iknbar=0; iknbar<Gk._size; iknbar++){
        for (size_t kkbar=0; kkbar<Gk._kArr_l.size(); kkbar++){
            Hubbard::K_1D kbar(Gk._kArr_l[kkbar],Gk._precomp_wn[iknbar]);
            upper_level += Gk(kbar._iwn,kbar._qx)(1,1) * gamma_oneD_spsp_full_middle(Gk,kbar,q) * Gk((kbar-q)._iwn,(kbar-q)._qx)(1,1);
        }
    }
    upper_level *= -1.0/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk);
    return upper_level;
}


/* Leaving the functions entering the full susceptibility. */

std::complex<double> Susceptibility::gamma_oneD_spsp(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D ktilde,Hubbard::K_1D kbar,Hubbard::K_1D q) const{ // q's contain bosonic Matsubara frequencies.
    std::complex<double> lower_level=0.0+0.0*im;
    for (int wttilde=0; wttilde<Gk._size; wttilde++){
        for (size_t qttilde=0; qttilde<Gk._kArr_l.size(); qttilde++){
            lower_level += Gk((ktilde+q)._iwn-Gk._precomp_qn[wttilde],(ktilde+q)._qx-Gk._kArr_l[qttilde])(0,0)*Gk((kbar+q)._iwn-Gk._precomp_qn[wttilde],(kbar+q)._qx-Gk._kArr_l[qttilde])(1,1);
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
                            Gk._precomp_wn[wtilde]+q._iwn,Gk._kArr_l[ktilde]+q._qx
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
            lower_level += Gk(wtilde+q._iwn-Gk._precomp_qn[wttilde],ktilde+q._qx-Gk._kArr_l[qttilde])(0,0)*Gk(wbar+q._iwn-Gk._precomp_qn[wttilde],kbar+q._qx-Gk._kArr_l[qttilde])(1,1);
        }
    }
    lower_level *= -1.0*Gk._u/(Gk._beta*Gk._Nk); /// Removed minus sign
    lower_level += 1.0;
    return Gk._u/lower_level;
}

std::tuple< std::complex<double>, std::complex<double> > Susceptibility::gamma_oneD_spsp_plotting(Hubbard::FunctorBuildGk& Gk,double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar,Hubbard::K_1D q) const{ // q's contain bosonic Matsubara frequencies.
    std::complex<double> lower_level=0.0+0.0*im, chi_bubble=0.0+0.0*im; // chi_bubble represents the lower bubble in the vertex function, not the total vertex function.
    for (int wttilde=0; wttilde<Gk._size; wttilde++){
        for (size_t qttilde=0; qttilde<Gk._kArr_l.size(); qttilde++){
            lower_level += Gk(wtilde+q._iwn-Gk._precomp_qn[wttilde],ktilde+q._qx-Gk._kArr_l[qttilde])(0,0)*Gk(wbar+q._iwn-Gk._precomp_qn[wttilde],kbar+q._qx-Gk._kArr_l[qttilde])(1,1);
        }
    }
    lower_level *= Gk._u/(Gk._beta*Gk._Nk); /// Removed minus sign
    chi_bubble = lower_level;
    lower_level *= -1.0;
    lower_level += 1.0;
    return std::make_tuple(Gk._u/lower_level,chi_bubble);
}

#ifdef PARALLEL

std::complex<double> Susceptibility::chispsp_long_expr(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const{
    std::ofstream output;
    std::string strOutput("ktilde_kbar_U"+std::to_string(Gk._u)+"_beta"+std::to_string(Gk._beta)+"_ndo"+std::to_string(Gk._ndo)+"_Nk"+std::to_string(Gk._Nk)+"_Nw"+std::to_string(Gk._size)+"__AA.dat");
    std::complex<double> upper_level=0.0+0.0*im; // total susceptibility
    arma::Mat< std::complex<double> > ktb(Gk._kArr_l.size(),Gk._kArr_l.size()); // Container for ktilde-kbar susceptibility!
    unsigned int it=0;
    int totSize=Gk._kArr_l.size()*Gk._kArr_l.size();
    std::cout << "totSize: " << totSize << " size kArr_l: " << Gk._kArr_l.size() << std::endl;
    ThreadFunctor::ThreadFunctor1D threadObj(upper_level,Gk,q,ktb);
    while (it<totSize){
        if (totSize % NUM_THREADS != 0){
            if ( totSize-it == (totSize % NUM_THREADS) ){
                int newl=totSize-it;
                std::vector<std::thread> tt(newl);
                for (int l=0; l<newl; l++){
                    int ltot=it+l; // Have to make sure spans over the whole array of k-space.
                    int lkt = static_cast<int>(floor(ltot/Gk._kArr_l.size()));
                    int lkb = (totSize % Gk._kArr_l.size());
                    std::thread t(threadObj,lkt,lkb);
                    tt[l]=std::move(t);
                    std::cout << "ho" << ltot << "\n";
                }
                threadObj.join_all(tt);
            }
            else{
                std::vector<std::thread> tt(NUM_THREADS);
                for (int l=0; l<NUM_THREADS; l++){
                    int ltot=it+l; // Have to make sure spans over the whole array of k-space.
                    int lkt = static_cast<int>(floor(ltot/Gk._kArr_l.size()));
                    int lkb = (ltot % Gk._kArr_l.size());
                    std::cout << "lkt: " << lkt << " lkb: " << lkb << "\n";
                    std::thread t(threadObj,lkt, lkb);
                    tt[l]=std::move(t);
                    std::cout << "hola" << ltot << "\n";
                }
                threadObj.join_all(tt);
            }
        }
        else{
            std::vector<std::thread> tt(NUM_THREADS);
            for (int l=0; l<NUM_THREADS; l++){
                int ltot = it+l; // Have to make sure spans over the whole array of k-space.
                int lkt = static_cast<int>(floor(ltot/Gk._kArr_l.size()));
                int lkb = (ltot % Gk._kArr_l.size());
                std::thread t(threadObj,lkt, lkb);
                tt[l]=std::move(t);
            }
            threadObj.join_all(tt);
        }
        it+=NUM_THREADS;
    }
    upper_level *= -1.0*(1.0/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk)/(Gk._beta*Gk._Nk));
    output.open(strOutput, std::ofstream::out | std::ofstream::app);
    for (int ktilde=0; ktilde<Gk._kArr_l.size(); ktilde++){ // Printing to file
        for (int kbar=0; kbar<Gk._kArr_l.size(); kbar++){
            output << ktb(ktilde,kbar) << " ";
        }
        output << "\n";
    }
    output.close();
    return upper_level;
}

#else

std::complex<double> Susceptibility::chispsp_long_expr(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const{
    std::ofstream output;
    std::string strOutput("ktilde_kbar_U"+std::to_string(Gk._u)+"_beta"+std::to_string(Gk._beta)+"_ndo"+std::to_string(Gk._ndo)+"_Nk"+std::to_string(Gk._Nk)+"_Nw"+std::to_string(Gk._size)+"__AA.dat");
    std::complex<double> upper_level=0.0+0.0*im;
    for (size_t ktilde=0; ktilde<Gk._kArr_l.size(); ktilde++){
        std::cout << "ktilde: " << ktilde << std::endl;
        for (size_t kbar=0; kbar<Gk._kArr_l.size(); kbar++){
            std::cout << "kbar: " << kbar << std::endl;
            std::complex<double> tmp_val_kbar=0.0+0.0*im;
            // #pragma omp parallel for shared(Gk,ktilde,kbar,q) reduction(+: upper_level, tmp_val_kbar)
            for (int wtilde=0; wtilde<Gk._size; wtilde++){
                for (int wbar=0; wbar<Gk._size; wbar++){
                    std::complex<double> tmp_val = Gk(
                            Gk._precomp_wn[wtilde],Gk._kArr_l[ktilde]
                            )(0,0)*Gk(
                            Gk._precomp_wn[wtilde]+q._iwn,Gk._kArr_l[ktilde]+q._qx
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

#endif

#ifdef INTEGRAL

std::complex<double> Susceptibility::chi0(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const{
    Integrals integralsObj;
    std::complex<double> xi=0.0+0.0*im;
    for (size_t ikn=0; ikn<Gk._size; ikn++){ // Integrating over BZ for each Matsubara frequency.
        std::function< double(double) > tmpFunctImag = [&](double k){
            return ( Gk( Gk._precomp_wn[ikn] + q._iwn, k + q._qx )(0,0) * Gk( Gk._precomp_wn[ikn], k )(0,0) ).imag();
        };
        std::function< double(double) > tmpFunctReal = [&](double k){
            return ( Gk( Gk._precomp_wn[ikn] + q._iwn, k + q._qx )(0,0) * Gk( Gk._precomp_wn[ikn], k )(0,0) ).real();
        };
        xi += ( integralsObj.integrate_simps(tmpFunctReal,-M_PI,M_PI,0.001) + im*integralsObj.integrate_simps(tmpFunctImag,-M_PI,M_PI,0.001) ); // 0.001 is enough!
    }
    xi *= 1.0/(Gk._beta*(2.0*M_PI)); /// Removed minus sign
    return xi;
}

#else

std::complex<double> Susceptibility::chi0(Hubbard::FunctorBuildGk& Gk,Hubbard::K_1D q) const{
    std::complex<double> xi=0.0+0.0*im;
    for (size_t ikn=0; ikn<Gk._size; ikn++){
        for (size_t k=0; k<Gk._kArr_l.size(); k++){
            xi += Gk( Gk._precomp_wn[ikn] + q._iwn, Gk._kArr_l[k] + q._qx )(0,0) * Gk( Gk._precomp_wn[ikn], Gk._kArr_l[k] )(0,0);
            //xi += Gk( Gk._precomp_wn[ikn] + q._iwn, Gk._kArr_l[k] + q._qx )(1,1) * Gk( Gk._precomp_wn[ikn], Gk._kArr_l[k] )(1,1); // Changed (0,0) for (1,1)
        }
    }
    xi *= -1.0/(Gk._beta*Gk._Nk); /// Included minus sign
    return xi;
}

#endif
        
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

std::tuple< std::complex<double>, std::complex<double> > Susceptibility::get_chi_1D(Hubbard::FunctorBuildGk& Gk, std::string filename_chi, std::string filename_chi0){
    /* Prints out in distinct files the imaginary values of the chi and chi0. */
    std::ofstream outputFileChi;
    std::ofstream outputFileChi0;
    std::complex<double> Suscep, Suscep0;
    std::complex<double> TVSUSusChi0(0.0,0.0);
    std::complex<double> TVSUSusChi(0.0,0.0);
    std::tuple< std::complex<double>, std::complex<double> > suss;
    for (int k=0; k<=Gk._Nk; k++){
        Suscep = chisp(Gk,Hubbard::K_1D(Gk._kArr_l[k],0.0+0.0*im)); // Bosonic Matsubara frequency!
        Suscep0 = chi0(Gk,Hubbard::K_1D(Gk._kArr_l[k],0.0+0.0*im)); // Bosonic Matsubara frequency!
        // if ((k==0) || (k==Gk._Nk)){ // only keeping susceptibility at k=pi and -pi
        //     TVSUSusChi+=Suscep;
        //     TVSUSusChi0+=Suscep0;
        // }
        if (k == (int)(Gk._Nk/4)){ // only keeping susceptibility at k=pi and -pi
            TVSUSusChi+=Suscep;
            TVSUSusChi0+=Suscep0;
        }
        /* Saving imaginary part of the susceptibility */
        outputFileChi.open(filename_chi, std::ofstream::out | std::ofstream::app);
        outputFileChi << Suscep.imag() << " ";
        outputFileChi.close();
        // cout << "susceptibility at k: " << kArr_l[k] << " is " << Suscep << endl;
        outputFileChi0.open(filename_chi0, std::ofstream::out | std::ofstream::app);
        outputFileChi0 << Suscep0.imag() << " ";
        outputFileChi0.close();
    }
    // TVSUSusChi0 *= 0.5;
    // TVSUSusChi *= 0.5;
    std::cout << "susceptibility chi: " << TVSUSusChi << " vs chi0: " << TVSUSusChi0 << std::endl;
    outputFileChi.open(filename_chi, std::ofstream::out | std::ofstream::app);
    outputFileChi << "\n";
    outputFileChi.close();
    outputFileChi0.open(filename_chi0, std::ofstream::out | std::ofstream::app);
    outputFileChi0 << "\n";
    outputFileChi0.close();
    suss = std::make_tuple(TVSUSusChi,TVSUSusChi0);
    return suss;
}