# -*- coding: utf-8 -*-
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from PIL import Image
from io import BytesIO

if __name__ == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    # ======================================================
    # color scheme
    cmap = mpl.colormaps['Blues']
    color1 = cmap(0.3)
    color2 = cmap(0.6)
    color3 = cmap(0.9)

    # ======================================================
    # growth parameter: theta_g = 1 + 0.05*t, t from 0 to 4.6
    growth = np.arange(1, 1.23125, 0.00125)

    # ======================================================
    # load wiring length data
    # eta = 1 (K1 = K2), varying stiffness
    #   K1_eff = K2_eff = 100 N/m^2
    #   K1_eff = K2_eff = 200 N/m^2
    #   K1_eff = K2_eff = 300 N/m^2
    with open('Job-K1-100-K2-100-total.npy', 'rb') as f:
        total_100 = np.load(f)
    with open('Job-K1-200-K2-200-total.npy', 'rb') as f:
        total_200 = np.load(f)
    with open('Job-K1-300-K2-300-total.npy', 'rb') as f:
        total_300 = np.load(f)

    # ======================================================
    # plot normalized total wiring length
    plt.figure()
    plt.plot(growth, total_100 / total_100[0], color=color1, linestyle='solid')
    plt.plot(growth, total_200 / total_200[0], color=color2, linestyle='solid')
    plt.plot(growth, total_300 / total_300[0], color=color3, linestyle='solid')

    plt.legend([
        r'$K_{1,\rm eff} = K_{2,\rm eff} = 100$ N/m$^2$',
        r'$K_{1,\rm eff} = K_{2,\rm eff} = 200$ N/m$^2$',
        r'$K_{1,\rm eff} = K_{2,\rm eff} = 300$ N/m$^2$',
    ])

    plt.xlabel(r"$\vartheta_{g}$ [-]")
    plt.ylabel(r"normalized length $\overline{L}(t)$ [-]")

    plt.rcParams['figure.dpi'] = 500
    plt.rcParams['savefig.dpi'] = 500
    plt.tight_layout()

    # save figure to TIFF
    png1 = BytesIO()
    plt.savefig(png1, format='png', dpi=500)
    png2 = Image.open(png1)
    png2.save('figure_total_wiring_length_keep_eta.tiff')
    png1.close()

    plt.savefig("figure_total_wiring_length_keep_eta.png", dpi=500)
    plt.show()
