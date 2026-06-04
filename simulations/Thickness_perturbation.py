# author: Xincheng Wang
# contributors: Shuolun Wang 
# unit system: mm-N-MPa unit system
# dimension: plane strain

# -*- coding: utf-8 -*-
from __future__ import print_function
import Python_subroutine_axon_tension
import numpy as np
import os

# preload all subroutines to the local python environment
execfile("Python_subroutine_axon_tension.py")

if __name__ == '__main__':

    # ======================================================
    # dimensions of the model 
    Length = 80.0 #mm
    Height = 40.0 #mm
    Cortex_thickness = 2.0 #mm
    Dimensions = [Length, Height, Cortex_thickness]

    # ======================================================
    # material properties 
    mu_subcortex = 100.0e-6 # MPa
    lame_subcortex = 9.3 * mu_subcortex # MPa
    stiffness_ratio = 3
    mu_cortex = stiffness_ratio * mu_subcortex  # MPa
    lame_cortex = 9.3 * mu_cortex # MPa
    growth_rate = 0.05
    Materials = [mu_subcortex, lame_subcortex, mu_cortex, lame_cortex, growth_rate]

    # ======================================================
    # step timing parameters, increased increment limits and reduced stabilization 
    # to ensure physical accuracy and minimize artificial energy dissipation.

    Totaltime = 4.6               # Total growth time
    Maxincnum = 100000            # Increased for better convergence handling
    Defaultstabilization = 1e-5   # Reduced dissipated energy fraction
    Defaultdampingratio = 0.01    # Reduced adaptive stabilization tolerance
    Mininc = 1e-5                 # Minimum time increment
    Incrementsize = 0.025         # Initial time increment
    Steppara = [Totaltime, Maxincnum, Defaultstabilization, Defaultdampingratio, Mininc, Incrementsize]

    # ======================================================
    # axon tract parabola shape parameters
    a_coeff = 1./30. #1/mm
    b_coeff = -7.0 #mm
    m_coeff = 8.0 #mm (shifted off-center)

    # only primary axon tract
    curve_num = 1

    # parameters that control the segments
    Geometric_length = 0.13 #mm
    Stretch_ratio = 2
    Axon_tract_property = [Geometric_length, Stretch_ratio]

    # influence radius
    InfluenceRadius = 1.0 #mm

    # ======================================================
    # axon tract effective stiffness per unit depth
    # K_eff, scaled by 1e-6 for mm-N-MPa system
    Stiffness_range_SI = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]  # N/m^2
    Stiffness_range = [k * 1e-6 for k in Stiffness_range_SI]        # MPa

    # ======================================================
    # thickness perturbation parameterss
    Waves = 20
    # Perturbation as fraction of cortex thickness (unitless)
    Perturbation_range = [0.0250, 0.0275, 0.0300, 0.0325, 0.0350,
                          0.0375, 0.0400, 0.0425, 0.0450, 0.0475, 0.0500]

    # ======================================================
    # naming of model parts
    PartName = 'Part-1'
    Step = 'Step-1'
    InstanceName = 'Part-1-1'
    UMAT = "./UMAT_axon_tension.f"

    for i in range(0, len(Perturbation_range)):
        for j in range(0, len(Stiffness_range)):

            ModelName = 'Model-perturbation%d-axonstiffness%d' % (i, j)
            JobName = 'Job-perturbation%d-axonstiffness%d' % (i, j)

            Perturbation_percentage = Perturbation_range[i]
            Stiffness_primary = Stiffness_range[j]

            Create_Bilayered_Rectangle(ModelName, PartName, Dimensions)
            Create_Material(ModelName, Materials)
            Create_Section(ModelName, PartName, Dimensions)
            Create_Assembly(ModelName, PartName, InstanceName)
            Create_Sets(ModelName, InstanceName, Dimensions)
            Create_Step(ModelName, Step, Steppara)
            Create_Contact(ModelName, Step)
            Create_Boundary_Conditions(ModelName, Step)
            Create_Mesh(ModelName, PartName, Dimensions)
            Create_Axon_Connection(a_coeff, b_coeff, m_coeff, curve_num, ModelName, InstanceName, Axon_tract_property, Stiffness_primary, InfluenceRadius, Dimensions)
            Create_mesh_node_sets(ModelName, InstanceName)
            Create_thickness_perturbation(ModelName, Length, Waves, Perturbation_percentage, Cortex_thickness)
            Modify_input(ModelName)
            Modify_input_for_initialize_growth_variable(ModelName)
            Create_Job(ModelName, JobName, UMAT)

    # ==============================================================
    # post-processing: compute mean squared displacement psi
    # Benchmark: maximum thickness perturbation, minimum axon stiffness
    ODB_Name_Bench = 'Job-perturbation%d-axonstiffness%d.odb' % (10, 0)
    NodesetName = 'TOPSURF_NODES'
    Variable = 'COORD'

    [x_coords_bench, y_coords_bench] = Post_processing_odbs(ODB_Name_Bench, NodesetName, Step, Variable)

    ODBtype = 1
    parameter1_range = Perturbation_range
    parameter2_range = Stiffness_range
    psi_array = Calculate_psi(ODBtype, NodesetName, Step, Variable, parameter1_range, parameter2_range, x_coords_bench, y_coords_bench)
    psi_array = np.flip(psi_array, axis=0)

    # ==============================================================
    # write data to csv for plotting
    outdir = "../results"
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    np.savetxt(os.path.join(outdir, "psi_array_thickness_perturbation.csv"), psi_array, delimiter=",")