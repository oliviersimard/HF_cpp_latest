import matplotlib.pyplot as plt
import numpy as np
import getopt, sys
import re

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rcParams["font.size"] = 16
plt.rcParams["text.latex.preamble"]=[r"\usepackage[charter]{mathdesign}\usepackage{amsmath}"]

try:
    opts, args = getopt.getopt(sys.argv[1:], "hfo", ["help", "file=", "option="])
except getopt.GetoptError as err:
    print(err)
    sys.exit(2)

for o, a in opts:
    print("opts: ", opts)
    if o in ("-h", "--help"):
        print("To plot the susceptibility, use --option s","\n")
        print("To plot plot n vs U, use --option n","\n")
        print("To plot the imaginary-time Green's function, use the program plotOnSpot.sh.")
        sys.exit()
    elif o in ("-f", "--file"):
        assert type(a) == str, "Wrong type of entry: must be string!"
        filename=a
    elif o in ("-o", "--option"):
        assert type(a) == str, "Wrong type of entry: must be string!"
        option=a
    else:
        assert False, "unhandled option."

data = np.genfromtxt(filename,dtype=float,delimiter=" ")

index_beta = filename.find("beta")
index_U = [m.start() for m in re.finditer("U",filename)][-1] 

## Have to change this depending on the file. Has been automated using REGEX.
beta_init = float(filename[index_beta:].split("_")[1]); beta_step = float(filename[index_beta:].split("_")[2]); beta_max = float(filename[index_beta:].split("_")[3])
u_init = float(filename[index_U:].split("_")[1]); u_step = float(filename[index_U:].split("_")[2]); u_max = float(filename[index_U:].split("_")[3]) 

print(data.shape)

len_u = data.shape[1]
len_beta = data.shape[0]
u_arr = np.arange(u_init,u_max+u_step,u_step,dtype=float)
print("u lengths: ", len_u, len(u_arr))
beta_arr = np.arange(beta_init,beta_max+beta_step,beta_step,dtype=int)
temperature_arr = 1.0/beta_arr
print("beta lengths: ", len_beta, len(beta_arr))

assert len(u_arr)==len_u and len(beta_arr)==len_beta, "Error in size of arrays!! Check data size." 

print("shape of data: ",data.shape)

fig, ax = plt.subplots()

color=iter(plt.cm.rainbow(np.linspace(0,6,len_u)))

if option == "s":
    u_vals = []
    for b in range(len_beta):
        max_sus = np.max(data[b,:])  # Find the maximum y value
        index_of_max_sus = np.where(data[b,:] == max_sus)
        u_vals_el = u_arr[index_of_max_sus]
        u_vals.append(u_vals_el)


if option == "s":
    for l in range(len_beta):
        ax.plot(u_arr,data[l,:],marker='s',markersize=3,color=next(color),label=r'$\beta={0:3.1f}$'.format(beta_arr[l]))

    ax.set_title(r'$\chi$ vs $U$ ($\beta \in$ {0:3.1f},{1:3.1f})'.format(beta_init,beta_max), fontsize=20,y=1.04,loc='center')
    ax.set_xlabel(r'$U$', fontsize=20)
    ax.set_ylabel(r'$\chi_{\text{sp}}$', fontsize=20)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
elif option == "n":
    for l in range(len_beta):
        ax.plot(u_arr,data[l,:],marker='s',markersize=3,color=next(color),label=r'$\beta={0:3.1f}$'.format(beta_arr[l]))

    ax.set_title(r'$n^{AA}_{\uparrow}$ vs $U$', fontsize=20,y=1.04,loc='center')
    ax.set_xlabel(r'$U$', fontsize=20)
    ax.set_ylabel(r'$n^{AA}_{\uparrow}$', fontsize=20)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
else:
    raise(ValueError("Check --help for the options."))

plt.show()

if option == "s":
    fig2, ax2 = plt.subplots()

    ax2.set_title(r"$T$ vs $U$ (1D)",fontsize=20,y=1.04,loc='center')
    ax2.set_xlabel(r"U",fontsize=20)
    ax2.set_ylabel(r"T",fontsize=20)
    ax2.plot(u_vals,temperature_arr,marker='o',markersize=5)

    plt.show()
