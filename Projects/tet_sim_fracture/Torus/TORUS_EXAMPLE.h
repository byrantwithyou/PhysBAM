//#####################################################################
// Copyright 2003, Ron Fedkiw, Neil Molino, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TORUS_EXAMPLE
//#####################################################################
#ifndef __TORUS_EXAMPLE__
#define __TORUS_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/GRAIN_BOUNDARIES.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
using namespace PhysBAM;

template <class T>
class TORUS_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    GRAIN_BOUNDARIES<T> *grain_boundaries;
    int fracture_bias_direction_smoothing_passes;
    
    TORUS_EXAMPLE()
        :initial_height(1),initial_orientation((T)3.14/2.0,VECTOR_3D<T>(1,0,0)),
        initial_velocity(0,-10,0),initial_angular_velocity(0,0,0),grain_boundaries(0),
        fracture_bias_direction_smoothing_passes(0)
    {
        restart_step_number=0;
        final_time=5;
        frame_rate=100;
        cfl_number=50;
        cg_iterations=250;
        youngs_modulus=5e5;poissons_ratio=(T).4;Rayleigh_coefficient=(T).02;
        strcpy(output_directory,"Torus/output");
        strcpy(input_file,"../../Public_Data/Tetrahedralized_Volumes/torus_12k.tet");
        gravity=10;
        push_out=true;

        // Fracture Parameters
        perform_fracture=true;
        perform_self_collision=true;
        fracture_threshold_1=(T)1e4;fracture_threshold_2=(T)5*fracture_threshold_1;fracture_threshold_3=(T)10*fracture_threshold_1;
        compressive_threshold_1=(T)-1e4;compressive_threshold_2=(T)-5e15;compressive_threshold_3=(T)-1e15;
        bias_stress=true;
       
        fracture_bias_stress_scaling=1;
        fracture_bias_magnitude=9e3;
        fracture_bias_propagation=9e6;

        // Initiation Points
        interpolation_fraction_threshhold=(T).03;
        number_of_fracture_initiation_points=4;
    }

    virtual ~TORUS_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    tetrahedralized_volume.particles.Update_Position();
    
    std::fstream input;input.open(input_file,std::ios::in|std::ios::binary);
    if(!input.is_open()){std::cout << "problem opening input torus file" << std::endl;exit(0);}
    tetrahedralized_volume.Read_Float(input);input.close();
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tetrhedra = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 

    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();
    tetrahedralized_volume.particles.Store_Mass();

    tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 
    tetrahedralized_volume.Set_Density(500);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);

    // Set Initiation Points
    tetrahedralized_volume.Update_Bounding_Box();
    initiation_point_positions=new ARRAY<VECTOR_3D<T> >(1,number_of_fracture_initiation_points);
    initiation_point_reference_seed_positions=new ARRAY<VECTOR_3D<T> >(1,number_of_fracture_initiation_points);
    T radius=0.3*(tetrahedralized_volume.bounding_box->xmax-tetrahedralized_volume.bounding_box->xmin);
    T center_y=(tetrahedralized_volume.bounding_box->ymax+tetrahedralized_volume.bounding_box->ymin)/(T)2.0;
    T center_x=(tetrahedralized_volume.bounding_box->xmax+tetrahedralized_volume.bounding_box->xmin)/(T)2.0;
    T center_z=(tetrahedralized_volume.bounding_box->zmax+tetrahedralized_volume.bounding_box->zmin)/(T)2.0;
    for(int i=0;i<number_of_fracture_initiation_points;i++){
        T theta=i*6.28/(number_of_fracture_initiation_points);
        (*initiation_point_reference_seed_positions)(i)=VECTOR_3D<T>(radius*sin(theta)+center_x,center_y,radius*cos(theta)+center_z);}
    initiation_point_radii=new ARRAY<T>(1,number_of_fracture_initiation_points);
    T eta=6.28/(number_of_fracture_initiation_points);
    ARRAY<T>::copy(sin(eta/(T)2.0)*radius,*initiation_point_radii);
    //ARRAY<T>::copy(.4,*initiation_point_radii);
}
//#####################################################################
// Function Initialize_Bias_Stress_Constants
//#####################################################################
virtual void Initialize_Bias_Stress_Constants(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,FRACTURE_TETRAHEDRALIZED_VOLUME<T>& rftv)

{
    TET_SIM_EMBEDDED_EXAMPLE<T>::Initialize_Bias_Stress_Constants(etv,rftv);

    rftv.fracture_bias_direction_coefficient=0.3;
    rftv.eigenvector_coefficient=0;
    rftv.fracture_bias_propagation_coefficient=0.7;
    
    delete grain_boundaries;
    grain_boundaries=new GRAIN_BOUNDARIES<T>(rftv,rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh);
    grain_boundaries->Set_Random_Fracture_Bias_Directions();
    grain_boundaries->Smooth_Fracture_Bias_Directions(fracture_bias_direction_smoothing_passes);

    int number_of_regions=number_of_fracture_initiation_points;

    // Initialize Grain Seeding points
    grain_boundaries->seed_node_index=new ARRAY<int>(1,0);
    grain_boundaries->seed_tet_index=new ARRAY<int>(1,0);

    // Assign Starting tets based on positions of initiation points.
    for(int r=0;r<number_of_regions;r++) {
        T min_distance=1e10;
        int closest_tet=0;
        for(int t=0;t<etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m;t++) {
            VECTOR_3D<T> centroid=etv.tetrahedralized_volume.Centroid(t);
            if((centroid-(*initiation_point_positions)(r)).Magnitude()<min_distance) {
                min_distance=(centroid-(*initiation_point_positions)(r)).Magnitude();
                closest_tet=t;}}
        grain_boundaries->seed_tet_index->Append(closest_tet);
    }

    // Seed Nodes
    for(int r=0;r<number_of_regions;r++) {
        T min_distance=1e10;
        int closest_node=0;
        for(int n=1;n<=etv.particles.array_collection->Size();n++){
            if((etv.particles.X(n)-(*initiation_point_positions)(r)).Magnitude()<min_distance){
                min_distance=(etv.particles.X(n)-(*initiation_point_positions)(r)).Magnitude();
                closest_node=n;}}
        grain_boundaries->seed_node_index->Append(closest_node);
    }

    //grain_boundaries->Fill_Node_Regions_With_Uniform_Vectors(number_of_regions);
    grain_boundaries->Fill_Regions_With_Uniform_Vectors(number_of_regions);
    grain_boundaries->Print_Number_In_Regions(number_of_regions);
    //for(int t=0;t<rftv.fracture_bias_direction.m;t++)
    //    std::cout << "rftv.fracture_bias_direction(" << t << ")=" << rftv.fracture_bias_direction(t) << std::endl;


    //grain_boundaries->Smooth_Fracture_Bias_Directions(3);
    //grain_boundaries->Print_Number_In_Regions(number_of_regions);

    /*
    for(int t=0;t<rftv.fracture_bias_direction.m;t++){
        if(grain_boundaries->Tetrahedron_Contains_Nodes_From_Different_Regions(t)){
            rftv.fracture_bias_stress_scaling(t)=1000;}}
    */
}
//#####################################################################
};
#endif
