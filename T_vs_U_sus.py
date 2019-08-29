import matplotlib.pyplot as plt
import numpy as np
import getopt, sys
import re
from scipy.optimize import curve_fit

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rcParams["font.size"] = 16
plt.rcParams["text.latex.preamble"]=[r"\usepackage[charter]{mathdesign}\usepackage{amsmath}"]

def fit_phase_diagram(u, a, b):
    return np.exp(-1.0/(a*u+b))

def linear_fit_function(u,a,b):
    return a*u+b

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

try:
    opts, args = getopt.getopt(sys.argv[1:], "hfo", ["help", "file=", "option="])
except getopt.GetoptError as err:
    print(err)
    sys.exit(2)

for o, a in opts:
    print("opts: ", opts)
    if o in ("-h", "--help"):
        print("To plot the susceptibility, use --option s","\n")
        print("To plot n vs U, use --option n","\n")
        print("To plot susceptibility vs k-points, use --option k","\n")
        print("To plot the imaginary-time Green's function, use the program plotOnSpot.sh.")
        sys.exit()
    elif o in ("-f", "--file"):
        assert type(a) == str, "Wrong type of entry: must be string!"
        filename=a
    elif o in ("-o", "--option"):
        assert type(a) == str, "Wrong type of entry: must be string!"
        option=a
    else:
        assert False, bcolors.FAIL+"unhandled option."+bcolors.ENDC

data = np.genfromtxt(filename,dtype=float,delimiter=" ")

index_beta = filename.find("beta")
index_Nk = filename.find("Nk")
index_Nomega = filename.find("Nomega")
index_U = [m.start() for m in re.finditer("U",filename)][-1]

end_of_file = filename[index_Nk:]
if index_Nk<index_Nomega:
    end_of_file = filename[index_Nk:].rstrip(".dat")
else:
    end_of_file = filename[index_Nomega:].rstrip(".dat")

imageDir="/Users/simardo/Documents/PhD/HF_cpp/Latex_docs/images/"
rough_phase_diagram=False

## Have to change this depending on the file. Has been automated using REGEX.
if option == "s" or option == "n":

    u_init = float(filename[index_U:].split("_")[1]); u_step = float(filename[index_U:].split("_")[2]); u_max = float(filename[index_U:].split("_")[3]) 
    print("u_init: ", u_init, "u_step: ", u_step, "u_max: ", u_max, "\n")
    beta_init = float(filename[index_beta:].split("_")[1]); beta_step = float(filename[index_beta:].split("_")[2]); beta_max = float(filename[index_beta:].split("_")[3])
    print("beta_init: ", beta_init, "beta_step: ", beta_step, "beta_max: ", beta_max, "\n")

    print(bcolors.OKBLUE+"data shape: "+bcolors.ENDC, data.shape)

    len_u = data.shape[1]
    len_beta = data.shape[0]
    u_arr = np.arange(u_init,u_max+u_step,u_step,dtype=float)
    print("u lengths: ", len_u, len(u_arr))
    beta_arr = np.arange(beta_init,beta_max+beta_step,beta_step,dtype=int)
    temperature_arr = 1./beta_arr
    print("beta lengths: ", len_beta, len(beta_arr))

    assert len(u_arr)==len_u and len(beta_arr)==len_beta, "Error in size of arrays!! Check data size." 

    print("shape of data: ",data.shape)

    fig, ax = plt.subplots()

    color=iter(plt.cm.rainbow(np.linspace(0,2,len_u)))

    if option == "s": ## This part is to be used if Im part chosen. Won't work otherwise. Mean_chi* files here.
        u_vals = []
        if rough_phase_diagram and "real" in filename:
            for b in range(len_beta):
                # max_sus = np.max(data[b,:])  # Find the maximum y value
                # index_of_max_sus = np.where(data[b,:] == max_sus)
                # u_vals_el = u_arr[index_of_max_sus]
                # u_vals.append(u_vals_el)
                try:
                    AF_u_cut = np.array([uval for uval in data[b,:] if uval < -0.001],dtype=float)
                    index_of_AF_sus = np.where(data[b,:] == AF_u_cut[0])
                    print("index: ", index_of_AF_sus)
                    u_vals_el = u_arr[index_of_AF_sus]
                    u_vals.append(u_vals_el)
                except IndexError as err:
                    print(err)
        elif not rough_phase_diagram and "real" in filename:
            for b in range(len_beta):
                u_vals.append(1.0/data[b,:][0]) # Simply using the fact that 1-U\chi_0=0 at the denominator.
        else:
            print(bcolors.WARNING+"Achtung: To print the phase diagram, you have to pass in the real part of the bubble susceptibility."+bcolors.ENDC)


    if option == "s":
        for l in range(len_beta):
            ax.plot(u_arr,data[l,:],marker='s',markersize=3,color=next(color),label=r'$\beta={0:3.1f}$'.format(beta_arr[l]))

        chi_val=""
        if "chio" in filename:
            chi_val="chio"
            ylabel=r'$\chi^0(\pi)$'
        else:
            chi_val="chi"
            ylabel=r'$\chi_{\text{sp,b}}(\pi)$'

        ImOrRe=""
        if "imag" in filename:
            ImOrRe=r"$\operatorname{Im}$"
        elif "real" in filename:
            ImOrRe=r"$\operatorname{Re}$"
        else:
            raise(ValueError("Can either be real or imaginary parts."))

        ax.grid(True)
        ax.set_title(r'RPA spin susceptibility ($\beta \in$ {0:3.1f},{1:3.1f})'.format(beta_init,beta_max), fontsize=20,y=1.04,loc='center')
        ax.set_xlabel(r'$U$', fontsize=20)
        ax.set_ylabel(ImOrRe+ylabel, fontsize=20)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.gcf().set_size_inches(12,9)
        plt.savefig(imageDir+chi_val+"_sp_vs_U_{0:2.1f}_{1:2.1f}_{2:2.1f}_beta_{3:2.1f}_{4:2.1f}_{5:2.1f}_".format(u_init,u_step,u_max,beta_init,beta_step,beta_max)+end_of_file+".pdf")
       
    elif option == "n":
        for l in range(len_beta):
            ax.plot(u_arr,data[l,:],marker='s',markersize=3,color=next(color),label=r'$\beta={0:3.1f}$'.format(beta_arr[l]))

        ax.set_title(r'$n^{AA}_{\uparrow}$ vs $U$', fontsize=20,y=1.04,loc='center')
        ax.set_xlabel(r'$U$', fontsize=20)
        ax.set_ylabel(r'$n^{AA}_{\uparrow}$', fontsize=20)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.gcf().set_size_inches(12,9)
        plt.savefig(imageDir+"n_AA_up_vs_U_{0:2.1f}_{1:2.1f}_{2:2.1f}_beta_{3:2.1f}_{4:2.1f}_{5:2.1f}_".format(u_init,u_step,u_max,beta_init,beta_step,beta_max)+end_of_file+".pdf")

    if option == "s":
        if "real" in filename:
            fig2, axs = plt.subplots(2,1,sharex=True)
            axs[0].grid(True)
            axs[1].grid(True)

            popt, pcov = curve_fit(linear_fit_function,u_vals,-1.0/np.log(temperature_arr),p0=[0.22,0.0],method="lm")

            plt.subplots_adjust(hspace=0.03)
            u_vals=np.asarray(u_vals)
            print("temperature_arr: ",type(temperature_arr)," ",type(u_vals))
            axs[0].set_title(r"$T$ vs $U$ (1D)",fontsize=20,y=1.04,loc='center')
            axs[0].set_ylabel(r"T",fontsize=20)
            axs[0].plot(u_vals,temperature_arr,marker='o',markersize=5)
            axs[0].plot(u_vals,fit_phase_diagram(u_vals,*popt),marker='o',markersize=5,color="red",label='fit exp: a={0:4.2f}, b={1:4.2f}'.format(*popt))
            axs[0].annotate(r"$T(U)\propto e^{-\frac{1}{aU+b}}$", xy=(1.2,0.15), xytext=(1.2,0.15))
            axs[0].legend()

            axs[1].plot(u_vals,-1.0/np.log(temperature_arr),marker='o',markersize=5)

            # aa=0.22; bb=-0.0
            # u_vals_fit=fit_phase_diagram(u_vals,aa,bb)
            axs[1].plot(u_vals,linear_fit_function(u_vals,*popt),marker='o',markersize=5,color="red",label='linear fit: a={0:4.2f}, b={1:4.2f}'.format(*popt))
            axs[1].set_ylabel(r"$-\frac{1}{\ln{T}}$",fontsize=20)
            axs[1].set_xlabel(r"U",fontsize=20)
            axs[1].annotate(r"$T(U)\propto aU+b$", xy=(1.2,0.5), xytext=(1.2,0.5))
            axs[1].legend()

            plt.gcf().set_size_inches(12,9)
            plt.savefig(imageDir+"T_vs_U_phase_diagram_U_{0:2.1f}_{1:2.1f}_{2:2.1f}_beta_{3:2.1f}_{4:2.1f}_{5:2.1f}_rough_diag_{6}_".format(u_init,u_step,u_max,beta_init,beta_step,beta_max,rough_phase_diagram)+end_of_file+".pdf")


elif option == "k":

    u_init = float(filename[index_U:].split("_")[1]); u_step = float(filename[index_U:].split("_")[2]); u_max = float(filename[index_U:].split("_")[3]) 
    print("u_init: ", u_init, "u_step: ", u_step, "u_max: ", u_max, "\n")
    beta = float(filename[index_beta:].split("_")[1])
    print("beta: ", beta, "\n")

    print("data shape: ", data.shape)

    len_k = data.shape[1]
    len_u = data.shape[0]
    u_arr = np.arange(u_init,u_max,u_step,dtype=float)
    print("u lengths: ", len_u, len(u_arr))

    assert len(u_arr)==len_u, "Error in size of arrays!! Check data size."

    chi_val=""
    if "Chi0" in filename:
        chi_val="chio"
        ylabel=r'$\chi^0(\mathbf{k})$'
    else:
        chi_val="chi"
        ylabel=r'$\chi_{\text{sp,b}}(\mathbf{k})$'

    im_or_re_val=""
    ImOrRe=""
    if "imag" in filename:
        im_or_re_val="IM_"
        ImOrRe=r"$\operatorname{Im}$"
    elif "real" in filename:
        im_or_re_val="RE_"
        ImOrRe=r"$\operatorname{Re}$"
    else:
        raise(ValueError("Can either be real or imaginary parts."))

    def format_func(value, tick_number):
        # find number of multiples of pi/4
        value = -np.pi + value*2.0*np.pi/(len_k-1) # Adapt to the length of k-space array.
        N = int(np.round(4 * value / np.pi))
        if N == 0:
            return "0"
        elif N == 1:
            return r"$\pi/4$"
        elif N == 2:
            return r"$\pi/2$"
        elif N == 4:
            return r"$\pi$"
        elif N % 4 > 0:
            if N % 2 == 0:
                return r"${0}\pi/2$".format(N // 2)
            else:
                return r"${0}\pi/4$".format(N)
        else:
            return r"${0}\pi$".format(N // 4)

    fig, ax = plt.subplots()

    color=iter(plt.cm.rainbow(np.linspace(0,1,len_u)))

    for l in range(len_u):
        ax.plot(data[l,:],marker='s',markersize=3,color=next(color),label=r'$U={0:3.2f}$'.format(u_arr[l]))

    ax.set_title(r'RPA spin susceptibility ($\beta={0:3}/t$)'.format(int(beta)),fontsize=20,y=1.04,loc='center')
    ax.set_xlabel(r'$\mathbf{k}$', fontsize=20)
    ax.set_ylabel(ImOrRe+ylabel, fontsize=20)
    ax.xaxis.set_tick_params(which='major',direction='inout',length=6)
    ax.yaxis.set_tick_params(which='major',direction='inout',length=6)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    # ax.xaxis.set_ticklabels([round(x,2) for x in np.linspace(-np.pi,np.pi,9)], color='black') ## Mapping the data onto Brillouin zone
    ax.set_xlim(left=0,right=len_k)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),ncol=2)

    plt.gcf().set_size_inches(15,10)
    plt.savefig(imageDir+"k_dependence_"+im_or_re_val+chi_val+"_U_{0:2.1f}_{1:2.1f}_{2:2.1f}_beta_{3:2.1f}_".format(u_init,u_step,u_max,beta)+end_of_file+".pdf")
    

else:
    raise(ValueError("Check --help for the options."))
