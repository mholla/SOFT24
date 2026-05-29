# author: Xincheng Wang
# contributors: Shuolun Wang 
# unit system: mm-N-MPa unit system
# dimension: plane strain

# -*- coding: utf-8 -*-
from __future__ import print_function
import Python_subroutine_axon_tension
import numpy as np

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
    # Note: beta = 5 here (increased from 3 to amplify wiring length change)
    mu_subcortex = 100.0e-6 # MPa
    lame_subcortex = 9.3 * mu_cortex # MPa
    stiffness_ratio = 3
    mu_cortex = stiffness_ratio * mu_subcortex  # MPa
    lame_cortex = 9.3 * mu_subcortex # MPa
    growth_rate = 0.05
    Materials = [mu_subcortex, lame_subcortex, mu_cortex, lame_cortex, growth_rate]

    # ======================================================
    # step timing parameters, increased increment limits and reduced stabilization 
    # to ensure physical accuracy and minimize artificial energy dissipation.

    Totaltime = 4.6               # Total growth time
    Maxincnum = 100000            # Increased for better convergence handling
    Defaultstabilization = 1e-5   # Reduced dissipated energy fraction
    Defaultdampingratio = 0.01    # Reduced adaptive stabilization tolerance
    Mininc = 0.001                # Minimum time increment
    Incrementsize = 0.025         # Initial time increment
    Steppara = [Totaltime, Maxincnum, Defaultstabilization, Defaultdampingratio, Mininc, Incrementsize]

    # number of frames for post-processing
    Increment_num = len(np.arange(0, Totaltime, Incrementsize)) + 1

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
    # Fixed primary: K1_eff = 100 N/m^2
    #   eta = 0.25 -> K2 = 400 N/m^2 (secondary stronger)
    #   eta = 1    -> K2 = 100 N/m^2 (equal)
    #   eta = 4    -> K2 = 25  N/m^2 (primary stronger)
    Stiffness_primary_SI = 100  # N/m^2
    Stiffness_primary = Stiffness_primary_SI * 1e-6  # MPa
    eta_range = [0.25, 1, 4]

    # ======================================================
    # naming of model parts
    PartName = 'Part-1'
    Step = 'Step-1'
    InstanceName = 'Part-1-1'
    UMAT = "./UMAT_axon_tension.f"
    NodeSetName_primary = "Attachment Points-1-Set-1"
    NodeSetName_secondary = "Attachment Points-1-Set-2"

    # ======================================================
    # preallocate data storage
    length_t_primary = np.zeros(Increment_num)
    length_t_secondary = np.zeros(Increment_num)
    total_length_t = np.zeros(Increment_num)

    for i in range(0, curve_num3):

        eta = eta_range[i]
        Stiffness_secondary = Stiffness_primary / eta

        # (e.g., 100e-6 N/mm is converted to 100 for file naming)
        ModelName = 'Model-K1-%d-K2-%d' % (int(round(Stiffness_primary * 1e6)), int(round(Stiffness_secondary * 1e6)))
        JobName = 'Job-K1-%d-K2-%d' % (int(round(Stiffness_primary * 1e6)), int(round(Stiffness_secondary * 1e6)))
        ODB_Name = JobName + '.odb'

        Create_Bilayered_Rectangle(ModelName, PartName, Dimensions)
        Create_Material(ModelName, Materials)
        Create_Section(ModelName, PartName, Dimensions)
        Create_Assembly(ModelName, PartName, InstanceName)
        Create_Sets(ModelName, InstanceName, Dimensions)
        Create_Step(ModelName, Step, Steppara)
        Create_Contact(ModelName, Step)
        Create_Boundary_Conditions(ModelName, Step)
        Create_Mesh(ModelName, PartName, Dimensions)
        Create_Axon_Connection(a_coeff, b_coeff, m_coeff, curve_num1, ModelName, InstanceName, Axon_tract_property, Stiffness_primary, InfluenceRadius, Dimensions)
        Create_Axon_Connection(a_coeff, b_coeff, m_coeff, curve_num2, ModelName, InstanceName, Axon_tract_property, Stiffness_secondary, InfluenceRadius, Dimensions)
        Create_Axon_Connection(a_coeff, b_coeff, m_coeff, curve_num3, ModelName, InstanceName, Axon_tract_property, Stiffness_secondary, InfluenceRadius, Dimensions)
        Create_mesh_node_sets(ModelName, InstanceName)
        Modify_input_for_wiring_nodeSet(ModelName)
        Create_Job(ModelName, JobName, UMAT)

        # ==============================================================
        # calculate wiring length at each frame
        for FrameNumber in range(0, Increment_num):
            length_t_primary[FrameNumber] = Calculate_wiring_length(ODB_Name, NodeSetName_primary, FrameNumber)
            length_t_secondary[FrameNumber] = Calculate_wiring_length(ODB_Name, NodeSetName_secondary, FrameNumber)
            # Weighted total: eta * L1 + 2 * L2
            total_length_t[FrameNumber] = eta * length_t_primary[FrameNumber] + 2 * length_t_secondary[FrameNumber]

        # ==============================================================
        # save wiring length data
        primary_wiring_length_name = '../results/' + JobName + '-primary.npy'
        secondary_wiring_length_name = '../results/' + JobName + '-secondary.npy'
        total_wiring_length_name = '../results/' + JobName + '-total.npy'

        with open(primary_wiring_length_name, 'wb') as f:
            np.save(f, length_t_primary)
        with open(secondary_wiring_length_name, 'wb') as f:
            np.save(f, length_t_secondary)
        with open(total_wiring_length_name, 'wb') as f:
            np.save(f, total_length_t)
