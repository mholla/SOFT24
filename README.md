# Axon Tension Paper Code

This repository contains the Abaqus simulation scripts, post-processing utilities, and plotting data used to reproduce the computational results for ["Axonal tension contributes to consistent fold placement" (Wang et al. 2024)](https://doi.org/10.1039/D4SM00129J). The code builds bilayered brain-tissue models with embedded axon-tract tension, evaluates geometry- and stiffness-dependent folding responses, and recreates the figures from the processed simulation data.

## Software Requirements

- Abaqus/CAE 2022 with a working Fortran compiler for the user material subroutine.
- Python environment for plotting.
- Python packages listed in `requirements.txt`.
- A LaTeX installation is recommended because the plotting scripts use Matplotlib with `usetex=True`.

Install the plotting dependencies with:

```bash
python -m pip install -r requirements.txt
```

The simulation scripts should be run through Abaqus, not a standard Python interpreter, because they rely on Abaqus/CAE APIs and compile the UMAT subroutine.

## Repository Schematic

```text
SOFT24/
|-- README.md
|-- requirements.txt
|-- simulations/
|   |-- UMAT_axon_tension.f
|   |-- Python_subroutine_axon_tension.py
|   |-- Parabola_geometry.py
|   |-- Thickness_perturbation.py
|   |-- Three_axon_tracts.py
|   |-- Wiring_length_change_eta.py
|   `-- Wiring_length_keep_eta.py
|-- results/
|   |-- data_critical-strain.csv
|   |-- buckling_times.npy
|   |-- psi_array_*.csv
|   `-- Job-*-total.npy
`-- plotting/
    |-- plot_critical_strain.py
    |-- plot_heatmap_*.py
    |-- plot_*_parametric.py
    |-- plot_wiring_length_*.py
    `-- figure_*.png
```

## Reproducing Simulations

Run the Abaqus jobs from the `simulations/` directory so the scripts can find `UMAT_axon_tension.f` and `Python_subroutine_axon_tension.py` through their existing relative paths.

```bash
cd simulations
abaqus cae noGUI=Parabola_geometry.py
abaqus cae noGUI=Thickness_perturbation.py
abaqus cae noGUI=Three_axon_tracts.py
abaqus cae noGUI=Wiring_length_change_eta.py
abaqus cae noGUI=Wiring_length_keep_eta.py
```

The scripts create Abaqus models/jobs, post-process output databases where needed, and write processed arrays to `../results/`.

| Script | Purpose | Main processed output |
| --- | --- | --- |
| `Parabola_geometry.py` | Sweeps axon-tract parabola geometry. | `results/psi_array_geometry.csv` |
| `Thickness_perturbation.py` | Sweeps cortical thickness perturbation and axon stiffness. | `results/psi_array_thickness_perturbation.csv` |
| `Three_axon_tracts.py` | Sweeps primary and secondary tract stiffness for a three-tract model. | `results/psi_array_three_curves.csv` |
| `Wiring_length_change_eta.py` | Computes total wiring length while varying stiffness ratio `eta`. | `results/Job-K1-100-K2-*-total.npy` |
| `Wiring_length_keep_eta.py` | Computes total wiring length with `eta = 1` while varying stiffness. | `results/Job-K1-*-K2-*-total.npy` |

## Reproducing Figures

Most figure scripts are designed to be run from the `plotting/` directory and read processed data from `../results/`.

```bash
cd plotting
python plot_critical_strain.py
python plot_heatmap_geometry.py
python plot_heatmap_thickness_perturbation.py
python plot_heatmap_three_curves.py
python plot_thickness_perturbation_parametric.py
python plot_three_curves_parametric.py
```

The wiring-length plotting scripts currently load `.npy` files from the current working directory. A direct way to run them with the committed data is:

```bash
cd results
python ../plotting/plot_wiring_length_change_eta.py
python ../plotting/plot_wiring_length_keep_eta.py
```

These commands write the generated wiring-length figures to the current directory. To save them inside `plotting/`, copy or link the required `Job-*-total.npy` files into `plotting/` before running the scripts there.

| Figure/data source | Plotting script |
| --- | --- |
| `results/data_critical-strain.csv`, `results/buckling_times.npy` | `plotting/plot_critical_strain.py` |
| `results/psi_array_geometry.csv` | `plotting/plot_heatmap_geometry.py` |
| `results/psi_array_thickness_perturbation.csv` | `plotting/plot_heatmap_thickness_perturbation.py` |
| `results/psi_array_perturbation.csv` | `plotting/plot_thickness_perturbation_parametric.py` |
| `results/psi_array_three_curves.csv` | `plotting/plot_heatmap_three_curves.py`, `plotting/plot_three_curves_parametric.py` |
| `results/Job-K1-100-K2-25-total.npy`, `results/Job-K1-100-K2-100-total.npy`, `results/Job-K1-100-K2-400-total.npy` | `plotting/plot_wiring_length_change_eta.py` |
| `results/Job-K1-100-K2-100-total.npy`, `results/Job-K1-200-K2-200-total.npy`, `results/Job-K1-300-K2-300-total.npy` | `plotting/plot_wiring_length_keep_eta.py` |

## Notes

- The finite-element model uses the mm-N-MPa unit system and plane strain assumptions, as noted in the simulation scripts.
- Abaqus simulation outputs such as `.odb`, `.inp`, `.dat`, `.msg`, `.sta`, and `.log` files are generated during reproduction and are not required for plotting the committed processed results.
- The committed `results/` arrays allow the plotting scripts to be rerun without repeating the full Abaqus simulations.
