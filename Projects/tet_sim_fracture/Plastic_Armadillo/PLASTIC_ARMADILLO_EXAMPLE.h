//#####################################################################
// Copyright 2003, Neil Molino, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLASTIC_ARMADILLO_EXAMPLE
//#####################################################################
#ifndef __PLASTIC_ARMADILLO_EXAMPLE__
#define __PLASTIC_ARMADILLO_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/GRAIN_BOUNDARIES.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
using namespace PhysBAM;

template <class T>
class PLASTIC_ARMADILLO_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int m_input,n_input,mn_input;
    T xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input;
    ARRAY<int> constrained_nodes_set1,constrained_nodes_set2,constrained_nodes_set3;
    GRAIN_BOUNDARIES<T> *grain_boundaries;
    int fracture_bias_direction_smoothing_passes;
    T arm_pull_velocity;
    
    PLASTIC_ARMADILLO_EXAMPLE()
        :initial_height((T)1.2),initial_orientation(1,VECTOR_3D<T>(0,1,0)),
        initial_velocity(0,(T)-sqrt((T)80),0),initial_angular_velocity(0,0,0),
        m_input(2),n_input(8),mn_input(8),
        xmin_input((T)-.04),xmax_input((T)-xmin_input),ymin_input((T)-.5),ymax_input((T)-ymin_input),zmin_input((T)-.5),zmax_input((T)-zmin_input),
        grain_boundaries(0),fracture_bias_direction_smoothing_passes(0),arm_pull_velocity(1.0)
    {
        restart_step_number=0;
        final_time=20;
        frame_rate=60;
        cfl_number=20;
        cg_iterations=250;
        youngs_modulus=(T)5e3;poissons_ratio=(T).4;Rayleigh_coefficient=(T).02;
        strcpy(output_directory,"PLASTIC_ARMADILLO/output");
        cg_tolerance=1e-4;
        
        // Fracture Parameters
        substeps_before_build_etv=15;
        perform_fracture=false;
        fracture_threshold_1=(T)1e3;fracture_threshold_2=(T)5e3;fracture_threshold_3=(T)1e4;
        compressive_threshold_1=(T)-1e12;compressive_threshold_2=(T)-5e12;compressive_threshold_3=(T)-1e12;
        interpolation_fraction_threshhold=(T).1;
        number_of_fracture_initiation_points=4;
        fracture_bias_propagation=5e10;

        // Initiation Points
        initiation_point_positions=new ARRAY<VECTOR_3D<T> >(1,number_of_fracture_initiation_points);
        initiation_point_radii=new ARRAY<T>(1,number_of_fracture_initiation_points);
        ARRAY<T>::copy(.2,*initiation_point_radii);

        // Collision Parameters
        perform_self_collision=false;
        push_out=false;
        perturb_amount_for_collision_freeness=(T)1e-5;
        collisions_repulsion_spring_constant_over_mass_times_length=22000;
        floor_friction=.5;

        // Plasticity Parameters
        use_plasticity=true;
        yield_ratio=30;plastic_clamp_ratio=.5;

        // Other params
        solids_parameters.use_constant_mass=false;
        gravity=1;
    }

    virtual ~PLASTIC_ARMADILLO_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    std::fstream input;char filename[256];
    tetrahedralized_volume.particles.Update_Position();
    sprintf(filename,"../../Public_Data/Tetrahedralized_Volumes/armadillo_4K.tet");input.open(filename,std::ios::in|std::ios::binary);assert(input.is_open());
    tetrahedralized_volume.Read_Float(input);input.close();
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tetrhedra = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 

    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();
    tetrahedralized_volume.particles.Store_Mass();

    tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();

    tetrahedralized_volume.Update_Bounding_Box();
    T scale=0.02;
    T velocity=0;T miny=1e10;

    // Adjust for orientation.
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());
    T width=tetrahedralized_volume.bounding_box->zmax-tetrahedralized_volume.bounding_box->zmin;

    // Angular velocity and orientation
    for(int p=1;p<=tetrahedralized_volume.particles.array_collection->Size();p++){
        tetrahedralized_volume.particles.V(p)=VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(p)-center);
        tetrahedralized_volume.particles.X(p)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(p)-center);}

    // Move to ground instead of dropping
    for(int p=1;p<=tetrahedralized_volume.particles.array_collection->Size();p++) if(tetrahedralized_volume.particles.X(p).y<miny) miny=tetrahedralized_volume.particles.X(p).y;
    for(int p=1;p<=tetrahedralized_volume.particles.array_collection->Size();p++){
        tetrahedralized_volume.particles.X(p)=scale*(tetrahedralized_volume.particles.X(p)-VECTOR_3D<T>(0,miny,0));
        tetrahedralized_volume.particles.V(p)+=VECTOR_3D<T>(0,-velocity,0);}

    // Recompute Center
    tetrahedralized_volume.Update_Bounding_Box();
    center=tetrahedralized_volume.bounding_box->Center();
    std::cout << "Center: "<< tetrahedralized_volume.bounding_box->Center().x << " " << tetrahedralized_volume.bounding_box->Center().y <<  " " <<  tetrahedralized_volume.bounding_box->Center().z << std::endl;    
    width=tetrahedralized_volume.bounding_box->zmax-tetrahedralized_volume.bounding_box->zmin;

    // Constrain Nodes
    for(int p=1;p<=tetrahedralized_volume.particles.array_collection->Size();p++){
        if(tetrahedralized_volume.particles.X(p).y>center.y){
            if(tetrahedralized_volume.particles.X(p).z>center.z+width*0.3) constrained_nodes_set1.Append(p);
            else if(tetrahedralized_volume.particles.X(p).z<center.z-width*0.3) constrained_nodes_set2.Append(p);
            else if(tetrahedralized_volume.particles.X(p).y>center.y+0.7*(center.y-miny)) constrained_nodes_set3.Append(p);
        }
    }

    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
}
//#####################################################################
// Function Initialize_Bias_Stress_Constants
//#####################################################################
virtual void Initialize_Bias_Stress_Constants(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& eta,
                                              FRACTURE_TETRAHEDRALIZED_VOLUME<T>& rftv)
{
    if(perform_fracture) {
        rftv.Set_Fracture_Bias_Stress_Scaling(1);
        rftv.Set_Fracture_Bias_Magnitude(0);
        rftv.Set_Fracture_Bias_Propagation(0);

        rftv.fracture_bias_direction_coefficient=1;
        rftv.eigenvector_coefficient=0;
        rftv.fracture_bias_propagation_coefficient=0;
        
        delete grain_boundaries;
        grain_boundaries=new GRAIN_BOUNDARIES<T>(rftv,rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh);
        grain_boundaries->Set_Random_Fracture_Bias_Directions();
        grain_boundaries->Smooth_Fracture_Bias_Directions(fracture_bias_direction_smoothing_passes);

        int number_of_regions=5;
        //grain_boundaries->Fill_Node_Regions_With_Uniform_Vectors(number_of_regions);
        grain_boundaries->Fill_Regions_With_Uniform_Vectors(number_of_regions);
        grain_boundaries->Print_Number_In_Regions(number_of_regions);
        for(int t=1;t<=rftv.fracture_bias_direction.m;t++)
            std::cout << "rftv.fracture_bias_direction(" << t << ")=" << rftv.fracture_bias_direction(t) << std::endl;

        grain_boundaries->Smooth_Fracture_Bias_Directions(9);
        //grain_boundaries->Print_Number_In_Regions(number_of_regions);

        /*
        for(int t=1;t<=rftv.fracture_bias_direction.m;t++){
            if(grain_boundaries->Tetrahedron_Contains_Nodes_From_Different_Regions(t)){
                rftv.fracture_bias_stress_scaling(t)=1000;}}
        */
    }
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
virtual void Set_External_Velocities(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=1;i<=constrained_nodes_set1.m;i++) V(constrained_nodes_set1(i))=VECTOR_3D<T>(0,0,arm_pull_velocity);
    for(int i=1;i<=constrained_nodes_set2.m;i++) V(constrained_nodes_set2(i))=VECTOR_3D<T>(0,0,-arm_pull_velocity);
    for(int i=1;i<=constrained_nodes_set3.m;i++) V(constrained_nodes_set3(i))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
virtual void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=1;i<=constrained_nodes_set1.m;i++) V(constrained_nodes_set1(i))=VECTOR_3D<T>(0,0,0);
    for(int i=1;i<=constrained_nodes_set2.m;i++) V(constrained_nodes_set2(i))=VECTOR_3D<T>(0,0,0);
    for(int i=1;i<=constrained_nodes_set3.m;i++) V(constrained_nodes_set3(i)).y=(T)0;
} 
};
#endif
