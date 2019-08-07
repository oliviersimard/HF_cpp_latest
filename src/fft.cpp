#include "fft.h"


void FFT::fft_t2w_G(Hubbard::FunctorBuildGk param, std::vector<double> &y, std::vector<std::complex<double> > &z){
  
  std::complex<double>* in=new std::complex<double> [param._size];
  std::complex<double>* out=new std::complex<double> [param._size];
  fftw_plan p;
  for(int j=0;j<param._size;j++){
    in[j]=y[j]*exp(im*(double)j*M_PI/(double)param._size);
  }
  p=fftw_plan_dft_1d(param._size, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  double delta = param._beta/(2.0*(double)param._size);
  for(int k=0;k<param._size;k++){
    std::complex<double> w = std::complex<double>(0.0,(2.0*(double)k+1.0)*M_PI/param._beta);
    std::cout << "out[" << k << "]: " << sin(w.imag()*delta) << std::endl;
    z[k]=1.0/(w) - im*sin(w.imag()*delta)*exp(-w*delta)/delta/w/w - out[k]*(2.0/delta/w/w)*sin(w.imag()*delta)*sin(w.imag()*delta);
  }
  delete [] in;
  delete [] out;
  fftw_destroy_plan(p);
}

void FFT::fft_t2w(Hubbard::FunctorBuildGk param, std::vector<double> &y, std::vector<std::complex<double> > &z)
/*----------------------------------------------------------------------------
    This subroutine computes the Fourier transformation

       z(k) = sum_{j=0}^{n} w_{n,j}*y(j)*exp(i*w_{k}*tau_{j})*dtau.
       (k = 0, 1,..., n-1)

       tau_{j} = j*dtau     (dtau = beta/n)

       w_{k} = (2*k-n+1)*pi/beta

       w_{n,j} = 1/2   for j=0, n
               = 1     for 1<=j<=n-1
----------------------------------------------------------------------------*/
{
  std::complex<double>* in=new std::complex<double> [param._size];
  std::complex<double>* out=new std::complex<double> [param._size];
  fftw_plan p;
  in[0]=0.5*y[0];
  for(int j=1;j<param._size;j++){
    in[j]=y[j]*exp(-im*(double)(param._size-1)*(double)j*M_PI/(double)param._size);
  }

  p=fftw_plan_dft_1d(param._size, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  for(int k=0;k<param._size;k++){
    //std::cout << "out[" << k << "]" << out[k] << std::endl;
    z[k]=(out[k]+0.5*y[param._size]*exp(-im*(double)(param._size-1)*M_PI))*(double)param._beta/(double)param._size;
  }
  delete [] in;
  delete [] out;
  fftw_destroy_plan(p);
  
}

void FFT::fft_t2w_notc(Hubbard::FunctorBuildGk param, std::vector<double> &y, std::vector<std::complex<double> > &z)
/*----------------------------------------------------------------------------
    This subroutine computes the Fourier transformation

       z(k) = sum_{j=0}^{n} w_{n,j}*y(j)*exp(i*w_{k}*tau_{j})*dtau.
       (k = 0, 1,..., n-1)

       tau_{j} = j*dtau     (dtau = beta/n)

       w_{k} = (2*k-n+1)*pi/beta

       w_{n,j} = 1/2   for j=0, n
               = 1     for 1<=j<=n-1
----------------------------------------------------------------------------*/
{
  std::complex<double>* in=new std::complex<double> [param._size];
  std::complex<double>* out=new std::complex<double> [param._size];
  fftw_plan p;
  in[0]=0.5*y[0];
  for(int j=1;j<param._size;j++){
    in[j]=y[j]*exp(im*(double)j*M_PI/(double)param._size);
  }

  p=fftw_plan_dft_1d(param._size, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  for(int k=0;k<param._size;k++){
    //std::cout << "out[" << k << "]" << out[k] << std::endl;
    z[k]=(out[k]+0.5*y[param._size]*exp(im*M_PI))*(double)param._beta/(double)param._size;
  }
  delete [] in;
  delete [] out;
  fftw_destroy_plan(p);
  
}


void FFT::fft_w2t(Hubbard::FunctorBuildGk param, std::vector<std::complex<double> > &z, std::vector<double> &y)
/*----------------------------------------------------------------------------
    This subroutine computes the Fourier transformation

    y(j) = 1/beta*sum_{k=0}^{n-1} z(k)*exp(-i*w_{k}*tau_{j}).
    (j = 0, 1,..., n)

    tau_{j} = j*dtau     (dtau = beta/n)

    w_{k} = (2*k-n+1)*pi/beta
----------------------------------------------------------------------------*/
{
  std::complex<double>* in=new std::complex<double> [param._size];
  std::complex<double>* out=new std::complex<double> [param._size];
  fftw_plan p;
  for(int k=0;k<param._size;k++){
    in[k]=z[k];
  }
  p=fftw_plan_dft_1d(param._size, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  for(int j=0;j<param._size;j++){
    y[j]=(out[j]*exp(im*(double)(param._size-1)*(double)j*M_PI/(double)param._size)).real()/param._beta;
  }
  y[param._size]=0.0;
  for(int k=0;k<param._size;k++){
    y[param._size]+=(z[k]*exp(im*(double)(param._size-1)*M_PI)).real()/param._beta;
  }
  delete [] in;
  delete [] out;
  fftw_destroy_plan(p);
}

void FFT::fft_w2t_notc(Hubbard::FunctorBuildGk param, std::vector<std::complex<double> > &z, std::vector<double> &y)
/*----------------------------------------------------------------------------
    This subroutine computes the Fourier transformation

    y(j) = 1/beta*sum_{k=0}^{N-1} z(k)*exp(-i*w_{k}*tau_{j}).
    (j = 0, 1,..., N)

    tau_{j} = j*dtau     (dtau = beta/N)

    w_{k} = (2*k+1)*pi/beta
----------------------------------------------------------------------------*/
{
  std::complex<double>* in=new std::complex<double> [param._size];
  std::complex<double>* out=new std::complex<double> [param._size];
  fftw_plan p;
  for(int k=0;k<param._size;k++){
    in[k]=z[k];
  }
  p=fftw_plan_dft_1d(param._size, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  for(int j=0;j<param._size;j++){
    y[j]=(out[j]*exp(-im*(double)j*M_PI/(double)param._size)).real()/param._beta;
  }
  y[param._size]=0.0;
  for(int k=0;k<param._size;k++){
    y[param._size]+=(z[k]*exp(-im*M_PI)).real()/param._beta;
  }
  for(int j=0; j<=param._size; j++){
    y[j] *= 2.0;
  }
  delete [] in;
  delete [] out;
  fftw_destroy_plan(p);
}

std::vector<double> FFT::get_gtau1D(Hubbard::FunctorBuildGk fftObj){
  std::vector< std::complex<double> > z(fftObj._size,std::complex<double>(0.0,0.0));
  std::vector<double> y(fftObj._size+1,0.0);
  
  for (int jj=0; jj<fftObj._size; jj++){
    std::complex<double> n_k(0.0,0.0);
    for (int kk=0; kk<=fftObj._Nk; kk++){
      if ( (kk==0) || (kk==fftObj._Nk) ){
        n_k += 0.5*fftObj.buildGkAA_1D(jj,fftObj._mu,fftObj._beta,fftObj._u,fftObj._ndo,fftObj._kArr[kk])(0,0);
      }
      else{
        n_k += fftObj.buildGkAA_1D(jj,fftObj._mu,fftObj._beta,fftObj._u,fftObj._ndo,fftObj._kArr[kk])(0,0); // This way the pointer to the Green's function in the class is not modified in-place.
      }
    }
    n_k /= fftObj._Nk;
    z[jj] = n_k - 1./fftObj.w(jj,0.0,fftObj._beta);
  }

  fft_w2t_notc(fftObj,z,y);

  for (int jj=0; jj<=fftObj._size; jj++)
    y[jj] -= 0.5;

  return y;
}

std::vector<double> FFT::get_gtau2D(Hubbard::FunctorBuildGk fftObj){
  std::vector< std::complex<double> > z(fftObj._size,std::complex<double>(0.0,0.0));
  std::vector<double> y(fftObj._size+1,0.0);
  
  for (int jj=0; jj<fftObj._size; jj++){
    std::complex<double> n_k(0.0,0.0);
    for (int kkx=0; kkx<=fftObj._Nk; kkx++){
      for (int kky=0; kky<=fftObj._Nk; kky++){
        if ( (kky==0) || (kky==fftObj._Nk) || (kkx==0) || (kkx==fftObj._Nk) ){
          n_k += 0.5*fftObj.buildGkAA_2D(jj,fftObj._mu,fftObj._beta,fftObj._u,fftObj._ndo,fftObj._kArr[kkx],fftObj._kArr[kky])(0,0);
          if ( ((kkx==0) || (kkx==fftObj._Nk)) && ((kky==0) || (kky==fftObj._Nk)) ){
            n_k += 0.25*fftObj.buildGkAA_2D(jj,fftObj._mu,fftObj._beta,fftObj._u,fftObj._ndo,fftObj._kArr[kkx],fftObj._kArr[kky])(0,0); // This way the pointer to the Green's function in the class is not modified in-place.
          }
        }
        else{
          n_k += fftObj.buildGkAA_2D(jj,fftObj._mu,fftObj._beta,fftObj._u,fftObj._ndo,fftObj._kArr[kkx],fftObj._kArr[kky])(0,0);
        }
      }
    }
    n_k /= (fftObj._Nk*fftObj._Nk);
    z[jj] = n_k - 1./fftObj.w(jj,0.0,fftObj._beta);
  }

  fft_w2t_notc(fftObj,z,y);

  for (int jj=0; jj<=fftObj._size; jj++)
    y[jj] -= 0.5;

  return y;
}

