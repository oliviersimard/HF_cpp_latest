// #include <boost/thread.hpp> // Define the Boost headers here, otherwise messes up std::complex<...> and it becomes ugly in vscode.
// #include <boost/chrono.hpp>
#include <thread>
#include <mutex>
#include <complex>
#include <chrono>
#include "green_utils.h"

#define NUM_THREADS 4

//#define PARALLEL

static std::mutex mutx;
extern arma::Mat< std::complex<double> > matGamma; // Matrices used in case parallel.
extern arma::Mat< std::complex<double> > matWeigths;
extern arma::Mat< std::complex<double> > matTotSus;

namespace ThreadFunctor{

template <typename Lockable>
class strict_lock{ // For external locking when parallelizing code! Not used for now.
public:
    typedef Lockable lockable_type;

    explicit strict_lock(lockable_type& obj) : obj_(obj) {
        obj.lock(); // locks on construction
    }
    strict_lock() = delete;
    strict_lock(const strict_lock&) = delete;
    strict_lock& operator=(const strict_lock&) = delete;

    ~strict_lock() { obj_.unlock(); } //  unlocks on destruction 

    bool owns_lock(const lockable_type* l) const noexcept{
      return l == &obj_;
    }
private:
    lockable_type& obj_;
};

class ThreadFunctor1D{
    public:
        ThreadFunctor1D(std::complex<double>& upper_level, Hubbard::FunctorBuildGk& Gk, Hubbard::K_1D& q, arma::Mat< std::complex<double> >& mat);
        ~ThreadFunctor1D(){};
        ThreadFunctor1D(ThreadFunctor1D&& rhs);
        void operator()(int ktilde, int kbar);
        std::complex<double> gamma_oneD_spsp(double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar);
        ThreadFunctor1D& operator=(const ThreadFunctor1D&)=delete;
        ThreadFunctor1D& operator=(ThreadFunctor1D&& other)=delete;
        ThreadFunctor1D(const ThreadFunctor1D&);
        inline void wait(int seconds){
            std::this_thread::sleep_for(std::chrono::seconds(seconds));
        }
        inline void lock() {
            mutx.lock();
        }
        inline void unlock() {
            mutx.unlock();
        }
        void join_all(std::vector<std::thread>& grp);
        

    private:
        std::complex<double>& _upper_level;
        Hubbard::FunctorBuildGk& _Gk;
        Hubbard::K_1D& _q;
        arma::Mat< std::complex<double> >& _ktb;
};

class ThreadWrapper{
    public:
        //ThreadWrapper(Hubbard::FunctorBuildGk& Gk, Hubbard::K_1D& q, arma::cx_dmat::iterator matPtr, arma::cx_dmat::iterator matWPtr);
        ThreadWrapper(Hubbard::FunctorBuildGk Gk,Hubbard::K_1D& q,double ndo_converged);
        ThreadWrapper(Hubbard::FunctorBuildGk Gk,Hubbard::K_2D& q,double ndo_converged);
        ~ThreadWrapper(){};
        void operator()(int ktilde, int kbar, double b); // 1D
        void operator()(int kbarx_m_tildex, int kbary_m_tildey); // 2D
        std::complex<double> gamma_oneD_spsp(double ktilde,std::complex<double> wtilde,double kbar,std::complex<double> wbar);
        std::complex<double> gamma_twoD_spsp(double kbarx_m_tildex,double kbary_m_tildey,std::complex<double> wtilde,std::complex<double> wbar);
        std::vector< std::complex<double> > buildGK1D(std::complex<double> ik, double k);
        std::vector< std::complex<double> > buildGK2D(std::complex<double> ik, double kx, double ky);
        void join_all(std::vector<std::thread>& grp);
    private:
        double _ndo_converged;
        Hubbard::FunctorBuildGk _Gk;
        Hubbard::K_1D& _q;
        //arma::cx_dmat::iterator _ktb;
        //arma::cx_dmat::iterator _ktbW;
};

} /* ThreadFunctor */