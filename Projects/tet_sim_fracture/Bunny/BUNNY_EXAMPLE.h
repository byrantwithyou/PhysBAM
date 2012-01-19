//#####################################################################
// Copyright 2003, Neil Molino, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BUNNY_EXAMPLE
//#####################################################################
#ifndef __BUNNY_EXAMPLE__
#define __BUNNY_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
using namespace PhysBAM;

template <class T>
class BUNNY_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    
    BUNNY_EXAMPLE()
                                :initial_height((T).21),initial_orientation(0,VECTOR_3D<T>(1,0,0)),
                                initial_velocity(0,-1,0),initial_angular_velocity(0,0,0)
    {
        restart_step_number=0;
        final_time=5;
        frame_rate=100;
        cfl_number=.5;
        cg_iterations=250;
        youngs_modulus=9000000;poissons_ratio=(T).4;Rayleigh_coefficient=(T).01;
        strcpy(output_directory,"Bunny/output");
        strcpy(input_file,"../../Public_Data/Tetrahedralized_Volumes/bunny.tet");
        gravity=10;
        
        // Fracture Parameters
        perform_fracture=true;
        fracture_threshold_1=(T)1e4;fracture_threshold_2=(T)5*fracture_threshold_1;fracture_threshold_3=(T)100*fracture_threshold_1;
        compressive_threshold_1=(T)-fracture_threshold_1;compressive_threshold_2=(T)-5*fracture_threshold_1;compressive_threshold_3=(T)-100*fracture_threshold_1;

        interpolation_fraction_threshhold=(T)0.2;
        number_of_fracture_initiation_points=1;

        perform_self_collision=false;
    }

    virtual ~BUNNY_EXAMPLE()
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
    std::cout << "bounding_box size=" << tetrahedralized_volume.bounding_box->Size() << std::endl;
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    T scale=10;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){tetrahedralized_volume.particles.X(i)*=scale;tetrahedralized_volume.particles.V(i)*=scale;}
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 
    tetrahedralized_volume.Set_Density(500);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
}
//#####################################################################
};
#endif
