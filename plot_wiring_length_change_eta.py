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
    # K1_eff = 100 N/m^2, varying eta
    #   eta = 4    -> K2_eff = 25  N/m^2
    #   eta = 1    -> K2_eff = 100 N/m^2
    #   eta = 0.25 -> K2_eff = 400 N/m^2
    with open('Job-K1-100-K2-100-total.npy', 'rb') as f:
        total_eta1 = np.load(f)
    with open('Job-K1-100-K2-25-total.npy', 'rb') as f:
        total_eta4 = np.load(f)
    with open('Job-K1-100-K2-400-total.npy', 'rb') as f:
        total_eta025 = np.load(f)

    # ======================================================
    # plot normalized total wiring length
    plt.figure()
    plt.plot(growth, total_eta4 / total_eta4[0], color=color1, linestyle='solid')
    plt.plot(growth, total_eta1 / total_eta1[0], color=color2, linestyle='solid')
    plt.plot(growth, total_eta025 / total_eta025[0], color=color3, linestyle='solid')

    plt.legend([
        r'$\eta = 4$, $K_{2,\rm eff} = 25$ N/m$^2$',
        r'$\eta = 1$, $K_{2,\rm eff} = 100$ N/m$^2$',
        r'$\eta = 1/4$, $K_{2,\rm eff} = 400$ N/m$^2$',
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
    png2.save('figure_total_wiring_length_change_eta.tiff')
    png1.close()

    plt.savefig("figure_total_wiring_length_change_eta.png", dpi=500)
    plt.show()
