import pandas as pd
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib as mpl
from PIL import Image
from io import BytesIO
import numpy as np

mpl.rcParams.update(mpl.rcParamsDefault)
plt.rcParams["font.family"] = "Times New Roman"

if __name__ == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    plt.figure()
    plt.axis([0, 30, 0, 0.5])
    plt.xlabel(r"stiffness ratio $\beta$ [-]")
    plt.ylabel(r"critical strain $\varepsilon_{\rm crit}$ [-]")

    # name colors
    cmap = matplotlib.colormaps['Blues']
    color_0 = cmap(0.2)
    color_10 = cmap(0.4)
    color_150 = cmap(0.6)
    color_300 = cmap(0.8)

    # read in and plot data from Holland 2018
    data = pd.read_csv("data_critical-strain.csv")
    betas_theory = data['beta'].tolist()
    eps_theory = data['strain'].tolist()

    plt.plot(betas_theory, eps_theory, color=color_0, linestyle="-")
    
    # read in buckling times
    buckling_times = np.load("buckling_times.npy", allow_pickle=True).item()

    betas = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
    dot_theta_g = 0.05 
    
    # calculate theta_g(t) = 1 + dot_theta_g * t
    theta_g_10 = 1.0 + dot_theta_g * buckling_times["t_10"]
    theta_g_150 = 1.0 + dot_theta_g * buckling_times["t_150"]
    theta_g_300 = 1.0 + dot_theta_g * buckling_times["t_300"]

    # calculate critical strain epsilon = 1 - 1/theta_g
    eps_10 =  1.0 - (1.0 / theta_g_10)
    eps_150 = 1.0 - (1.0 / theta_g_150)
    eps_300 = 1.0 - (1.0 / theta_g_300)

    plt.plot(betas, eps_10, color=color_10, marker='o', linestyle='')
    plt.plot(betas, eps_150, color=color_150, marker='o', linestyle='')
    plt.plot(betas, eps_300, color=color_300, marker='o', linestyle='')
    
    K_SI = [10, 150, 300]
    plt.legend([
        'theoretical',
        r'$K = %d\ \mathrm{N/m}^{2}$' % K_SI[0],
        r'$K = %d\ \mathrm{N/m}^{2}$' % K_SI[1],
        r'$K = %d\ \mathrm{N/m}^{2}$' % K_SI[2],
    ])

    plt.rcParams['figure.dpi'] = 500
    plt.rcParams['savefig.dpi'] = 500


    # save figure to TIFF if needed
    png1 = BytesIO()
    plt.savefig(png1, format='png')
    png2 = Image.open(png1)
    png2.save('figure_critical-strain.tiff')
    png1.close()

    plt.savefig("figure_critical-strain.png")
    plt.show()