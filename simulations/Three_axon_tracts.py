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
    # mu_cortex = 100 Pa = 100e-6 MPa
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

    Totaltime = 4.0               # Total growth time, total growth limited to theta_g = 1 + 0.05*4.0 = 1.2
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
    m_coeff = 0.0 #mm (center)

    # axon tract numbers
    curve_num1 = 1 # primary (center)
    curve_num2 = 2 # secondary (right)
    curve_num3 = 3 # secondary (left)

    # parameters that control the segments
    Geometric_length = 0.13 #mm
    Stretch_ratio = 2
    Axon_tract_property = [Geometric_length, Stretch_ratio]

    # influence radius
    InfluenceRadius = 1.0 #mm

    # ======================================================
    # axon tract effective stiffness per unit depth
    # K_eff range: 10 to 1000 N/m^2 (SI), scaled by 1e-6 for mm-N-MPa
    num = 10
    min_stiffness_SI = 10    # N/m^2
    max_stiffness_SI = 1000  # N/m^2
    Stiffness_range_SI = np.linspace(min_stiffness_SI, max_stiffness_SI, num)
    Stiffness_range = Stiffness_range_SI * 1e-6  # MPa (model units)

    # ======================================================
    # naming of model parts
    PartName = 'Part-1'
    Step = 'Step-1'
    InstanceName = 'Part-1-1'
    UMAT = "./UMAT_axon_tension.f"

    for i in range(0, len(Stiffness_range)):
        for j in range(0, len(Stiffness_range)):

            ModelName = 'Model-primary%d-secondary%d' % (i, j)
            JobName = 'Job-primary%d-secondary%d' % (i, j)

            Stiffness_primary = Stiffness_range[i]
            Stiffness_secondary = Stiffness_range[j]

            Create_Bilayered_Rectangle(ModelName, PartName, Dimensions)
            Create_Material(ModelName, Materials)
            Create_Section(ModelName, PartName, Dimensions)
            Create_Assembly(ModelName, PartName, InstanceName)
            Create_Sets(ModelName, InstanceName, Dimensions)
            Create_Step(ModelName, Step, Steppara)
            Create_Contact(ModelName, Step)
            Create_Boundary_Conditions(ModelName, Step)
            Create_Mesh(ModelName, PartName, Dimensions)

            # Primary axon tract (center)
            Create_Axon_Connection(a_coeff, b_coeff, m_coeff, curve_num1, ModelName, InstanceName, Axon_tract_property, Stiffness_primary, InfluenceRadius, Dimensions)
            # Secondary axon tract (right)
            Create_Axon_Connection(a_coeff, b_coeff, m_coeff, curve_num2, ModelName, InstanceName, Axon_tract_property, Stiffness_secondary, InfluenceRadius, Dimensions)
            # Secondary axon tract (left)
            Create_Axon_Connection(a_coeff, b_coeff, m_coeff, curve_num3, ModelName, InstanceName, Axon_tract_property, Stiffness_secondary, InfluenceRadius, Dimensions)

            Create_mesh_node_sets(ModelName, InstanceName)
            Modify_input(ModelName)
            Modify_input_for_initialize_growth_variable(ModelName)
            Create_Job(ModelName, JobName, UMAT)

    # ==============================================================
    # post-processing: compute mean squared displacement psi
    # Benchmark: max primary stiffness (K1=1000 N/m^2), min secondary (K2=10 N/m^2)
    ODB_Name_Bench = 'Job-primary%d-secondary%d.odb' % (9, 0)
    NodesetName = 'TOPSURF_NODES'
    Variable = 'COORD'

    [x_coords_bench, y_coords_bench] = Post_processing_odbs(ODB_Name_Bench, NodesetName, Step, Variable)

    ODBtype = 2
    parameter1_range = Stiffness_range
    parameter2_range = Stiffness_range
    psi_array = Calculate_psi(ODBtype, NodesetName, Step, Variable, parameter1_range, parameter2_range, x_coords_bench, y_coords_bench)
    psi_array = np.flip(psi_array, axis=0)

    # ==============================================================
    # write data to csv for plotting
    outdir = "../results"
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    np.savetxt(os.path.join(outdir, "psi_array_three_curves.csv"), psi_array, delimiter=",")