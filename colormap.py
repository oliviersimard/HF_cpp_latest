import numpy as np
import getopt, sys
import matplotlib.pyplot as plt
import os

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

try:
    opts, args = getopt.getopt(sys.argv[1:], "hf:v", ["help", "file="])
except getopt.GetoptError as err:
    print(err) 
    sys.exit(2)

output = None
verbose = False
for o, a in opts:
    if o == "-v":
        verbose = True
    elif o in ("-h", "--help"):
        sys.exit()
    elif o in ("-f", "--file"):
        assert type(a)==str, "Wrong type of input command line!"
        filename = a
    else:
        assert False, "unhandled option"

getSizeSquare = int(os.popen('head -n1 {0} | grep -o " " | wc -l'.format(filename)).read())
mesh = np.empty(shape=(getSizeSquare,getSizeSquare))
with open(filename) as f:
    for k1,line in enumerate(f):
        if k1<getSizeSquare:
            for k2,el in enumerate(line.split(" ")):
                if "(" in el.strip("\n"):
                    print("el: ", el)
                    tup = tuple(float(i) for i in el.strip('()').split(','))
                    mesh[k1,k2] = tup[1] ## Imaginary part of the susceptibility
        else:
            break

assert mesh.shape[0]==mesh.shape[1], "The data doesn't have \"square\" dimension"
max_val = mesh.shape[0]-1
def main():

    fig, ax = plt.subplots(1, 1, figsize=(9, 9))
    im = ax.imshow(mesh, aspect="auto", origin='lower', cmap=plt.get_cmap('viridis'))
    ax.set_title(r"$\operatorname{Im}\chi_{\text{sp}}$ (ladder diagrams)")
    ax.set_xlabel(r"$\tilde{k}$", fontsize=25, labelpad=10.0)
    ax.set_ylabel(r"$\bar{k}$", fontsize=25, labelpad=10.0)
    ax.set_xticks(np.linspace(0,max_val,10),minor=False)
    ax.set_yticks(np.linspace(0,max_val,10),minor=False)
    ax.xaxis.set_tick_params(which='major',direction='inout',length=6)
    ax.yaxis.set_tick_params(which='major',direction='inout',length=6)
    ax.xaxis.set_ticklabels([round(x,2) for x in np.linspace(-np.pi,np.pi,10)], color='black', fontsize=10) ## Mapping the data onto Brillouin zone
    ax.yaxis.set_ticklabels([round(x,2) for x in np.linspace(-np.pi,np.pi,10)], color='black', fontsize=10)
    ax.set_xlim(left=0,right=max_val)
    ax.set_ylim(bottom=0,top=max_val)
    fig.colorbar(im, ax=ax)
    plt.show()

    return None

if __name__ == "__main__":
    main()
