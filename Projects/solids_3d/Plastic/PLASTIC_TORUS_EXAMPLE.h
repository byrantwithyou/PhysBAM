//#####################################################################
// Copyright 2002, 2003, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLASTIC_TORUS_EXAMPLE
//#####################################################################
#ifndef __PLASTIC_TORUS_EXAMPLE__
#define __PLASTIC_TORUS_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "TET_SIM_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class PLASTIC_TORUS_EXAMPLE:public TET_SIM_EXAMPLE<T>
{
public:
    bool use_plasticity;
    DIAGONALIZED_PLASTICITY_MODEL_3D<T>* plasticity_model;
    T yield_ratio,plastic_clamp_ratio;
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    char input_file[256];

    PLASTIC_TORUS_EXAMPLE()
        :use_plasticity(true),plasticity_model(0),yield_ratio(2),plastic_clamp_ratio(2),
        initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(5,0,0)
    {
        final_time=5;
        restart_step_number=0;   
        cfl_number=(T)10;
        cg_tolerance=(T)1e-2;
        use_masses_and_springs=false;use_altitude_springs=false;
        use_diagonalized_fvm=true;
        use_linear_elasticity=false;use_neo_hookean=false;
        youngs_modulus=200000;poissons_ratio=(T).45;Rayleigh_coefficient=(T).01;
        use_plasticity=true;
        yield_ratio=1.5;
        plastic_clamp_ratio=4;
        strcpy(output_directory,"Torus/output");
        strcpy(input_file,"Torus/adaptive_torus_float.tet");
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=true;
        solids_parameters.collide_with_interior=false;
    }

    ~PLASTIC_TORUS_EXAMPLE()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    HEAVY_PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;

    std::fstream input;input.open(input_file,std::ios::in|std::ios::binary);tetrahedralized_volume.Read_Float(input);input.close();
    std::cout << "total vertices = " << particles.Size() << std::endl;std::cout << "total tetrhedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();tetrahedralized_volume.particles.Store_Mass(); // add back with the proper sizes

    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<particles.array_size;i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}
}
//#####################################################################
// Function Initialize_Diagonalized_Finite_Volume_Model
//#####################################################################
void Initialize_Diagonalized_Finite_Volume_Model(DEFORMABLE_OBJECT<T,VECTOR_3D<T> >& deformable_object,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    TET_SIM_EXAMPLE<T>::Initialize_Diagonalized_Finite_Volume_Model(deformable_object,tetrahedralized_volume);
    if(use_plasticity){
        plasticity_model=new DIAGONALIZED_PLASTICITY_MODEL_3D<T>(tetrahedralized_volume,yield_ratio,plastic_clamp_ratio);
        diagonalized_fvm->plasticity_model=plasticity_model;}
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
void Initialize_Rigid_Bodies(ARRAY<RIGID_BODY<TV>*>& rigid_bodies)
{  
    std::fstream input;char filename[256];
    
    // plane
    rigid_bodies.Resize(rigid_bodies.m+1);rigid_bodies(rigid_bodies.m)=new RIGID_BODY<TV>;
    // triangulated surface
    int index=triangulated_surface_list.Add_Triangulated_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.tri");input.open(filename,std::ios::in|std::ios::binary);
    triangulated_surface_list.triangulated_surface(index)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->Initialize_Triangulated_Surface(*triangulated_surface_list.triangulated_surface(index));
    // implicit surface
    index=implicit_surface_list.Add_Levelset_Implicit_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.phi");input.open(filename,std::ios::in|std::ios::binary);
    implicit_surface_list.implicit_surface(index)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->Initialize_Implicit_Surface(*implicit_surface_list.implicit_surface(index));
    // rigid body
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.rgd");input.open(filename,std::ios::in|std::ios::binary);
    rigid_bodies(rigid_bodies.m)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->coefficient_of_friction=(T).3;
}
//#####################################################################
};
}
#endif
