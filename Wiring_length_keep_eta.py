import Python_subroutine_axon_tension
import numpy as np

# preload all subroutines to the local python environment
execfile("Python_subroutine_axon_tension.py")

if __name__ == '__main__':

    # ======================================================
    # dimensions of the model
    Length = 0.080 #m
    Height = 0.040 #m
    Cortex_thickness = 0.002 #m
    Dimensions = [Length,Height,Cortex_thickness]

    # ======================================================
    # material properties; note the stiffness ratio is changed to 5 here
    mu_cortex = 100.0 #Pa
    lame_cortex= 9.3*mu_cortex #Pa
    stiffness_ratio = 5 
    mu_subcortex = stiffness_ratio*mu_cortex #Pa
    lame_subcortex = 9.3*mu_subcortex #Pa
    growth_rate = 0.05 
    Materials = [mu_cortex,lame_cortex,mu_subcortex,lame_subcortex,growth_rate]

    # ======================================================
    # step timing parameters
    Totaltime = 4.6 
    Maxincnum = 1000 
    Defaultstabilization = 0.0002 
    Defaultdampingratio = 0.05 
    Mininc = 0.001 
    Incrementsize = 0.025
    Steppara = [Totaltime,Maxincnum,Defaultstabilization,Defaultdampingratio,Mininc,Incrementsize]
    # specify parameters for post-processing 
    Increment_num = len(np.arange(0, Totaltime, Incrementsize)) + 1

    # ======================================================
    # axon tract parabola shape parameters
    a_coeff = 1./30.*1000 #/m
    b_coeff = -7./1000 #m
    m_coeff = 0 #m

    # axon tract numbers
    curve_num1 = 1
    curve_num2 = 2
    curve_num3 = 3

    # parameters that control the segments
    Geometric_length = 0.00013 #m
    Stretch_ratio = 2
    InfluenceRadius = 0.001 #m
    Axon_tract_property = [Geometric_length,Stretch_ratio]

    # ======================================================
    # vary axon tract stiffness
    Stiffness_range = [100,200,300]
    eta = 1

    # ======================================================
    # naming of model parts
    PartName = 'Part-1'
    Step = 'Step-1'
    InstanceName = 'Part-1-1'
    UMAT = "./UMAT_axon_tension.f"
    NodeSetName_primary = "Attachment Points-1-Set-1"
    NodeSetName_secondary = "Attachment Points-1-Set-2"

    # ======================================================
    # preallocate data storage: wiring length change of the primary tract, the secondary tract, and the weighted total length
    length_t_primary = np.zeros(Increment_num)
    length_t_secondary = np.zeros(Increment_num)
    total_length_t = np.zeros(Increment_num)

    for i in range(0,curve_num3):

        Stiffness_primary = Stiffness_range[i]
        Stiffness_secondary = Stiffness_primary/eta 

        ModelName =  'Model-K1-%d-K2-%d' %(Stiffness_primary,Stiffness_secondary)
        JobName = 'Job-K1-%d-K2-%d' %(Stiffness_primary,Stiffness_secondary)
        ODB_Name = JobName +'.odb'

        Create_Bilayered_Rectangle(ModelName, PartName, Dimensions)
        Create_Material(ModelName,Materials)
        Create_Section(ModelName, PartName,Dimensions)
        Create_Assembly(ModelName,InstanceName)
        Create_Sets(ModelName,InstanceName,Dimensions)
        Create_Step(ModelName, Step, Steppara)
        Create_Contact(ModelName, Step)
        Create_Boundary_Conditions(ModelName, Step)
        Create_Mesh(ModelName, PartName, Dimensions)
        Create_Axon_Connection(a_coeff,b_coeff,m_coeff,curve_num1,ModelName,InstanceName,Axon_tract_property,Stiffness_primary,InfluenceRadius,Dimensions)
        Create_Axon_Connection(a_coeff,b_coeff,m_coeff,curve_num2,ModelName,InstanceName,Axon_tract_property,Stiffness_secondary,InfluenceRadius,Dimensions)
        Create_Axon_Connection(a_coeff,b_coeff,m_coeff,curve_num3,ModelName,InstanceName,Axon_tract_property,Stiffness_secondary,InfluenceRadius,Dimensions)
        Create_mesh_node_sets(ModelName,InstanceName)
        Modify_input_for_wiring_nodeSet(ModelName)
        Create_Job(ModelName, JobName, UMAT)

        # ==============================================================
        # calculate the wiring length change with respect to simulation time
        for FrameNumber in range(0,Increment_num):
            length_t_primary[FrameNumber]  = Calculate_wiring_length(ODB_Name,NodeSetName_primary,FrameNumber)
            length_t_secondary[FrameNumber] = Calculate_wiring_length(ODB_Name,NodeSetName_secondary,FrameNumber)
            total_length_t[FrameNumber] = length_t_primary[FrameNumber] + length_t_secondary[FrameNumber]*2/eta 
    
        # ==============================================================
        # specify variable names for post-processing 
        primary_wiring_length_name = JobName + '-primary.npy'
        secondary_wiring_length_name = JobName + '-secondary.npy'
        total_wiring_length_name = JobName + '-total.npy'

        # ==============================================================
        # write data to the .npy files for plotting
        with open(primary_wiring_length_name, 'wb') as f:
            np.save(f, length_t_primary)
        with open(secondary_wiring_length_name, 'wb') as f:
            np.save(f, length_t_secondary)
        with open(total_wiring_length_name, 'wb') as f:

            np.save(f, total_length_t)
