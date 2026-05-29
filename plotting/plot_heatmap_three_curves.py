# -*- coding: utf-8 -*-
import pandas as pd
from matplotlib import rc
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image
from io import BytesIO
import matplotlib as mpl
import numpy as np

mpl.rcParams.update(mpl.rcParamsDefault)
plt.rcParams["font.family"] = "Times New Roman"


def heat_map(csv_file):

    # read in data
    data = pd.read_csv(csv_file, header=None)

    # stiffness range: 10 to 1000 N/m^2 (10 equally spaced values)
    stiffness_ticks = np.linspace(10, 1000, 10).astype(int).tolist()

    # specify primary (K1) and secondary (K2) axon tract stiffness [N/m^2]
    data.columns = stiffness_ticks
    data.index = stiffness_ticks[::-1]

    # plot heatmap
    g = sns.heatmap(data, cmap='Blues', square=True,
                    linewidths=1, linecolor='white',
                    cbar_kws={'label': r'mean squared displacement $\psi$ [mm]',
                              'location': 'bottom', 'shrink': 0.51})

    # axes
    g.tick_params(left=False, bottom=False)
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    g.set_xticklabels(g.get_xticklabels(), rotation=0)
    g.set_xlabel(r"secondary axon tract stiffness $K_{2,\rm eff}$ [N/m$^2$]")
    g.set_ylabel(r"primary axon tract stiffness $K_{1,\rm eff}$ [N/m$^2$]")

    plt.tight_layout()
    plt.rcParams['figure.dpi'] = 500
    plt.rcParams['savefig.dpi'] = 500

    # save figure to TIFF
    png1 = BytesIO()
    plt.savefig(png1, format='png')
    png2 = Image.open(png1)
    png2.save('figure_heatmap_three_curves.tiff')
    png1.close()

    plt.savefig("figure_heatmap_three_curves.png")
    plt.show()


if __name__ == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    heat_map("../results/psi_array_three_curves.csv")
