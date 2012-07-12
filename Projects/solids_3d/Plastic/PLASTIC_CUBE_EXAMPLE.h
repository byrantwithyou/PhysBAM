//#####################################################################
// Copyright 2002, 2003, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLASTIC_CUBE_EXAMPLE
//#####################################################################
#ifndef __PLASTIC_CUBE_EXAMPLE__
#define __PLASTIC_CUBE_EXAMPLE__

#include "TET_SIM_EXAMPLE.h"
#include <Forces_And_Torques/DIAGONALIZED_PLASTICITY_CONTROL_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_SIMPLE_PLASTICITY_3D.h>
namespace PhysBAM{

template<class T>
class PLASTIC_CUBE_EXAMPLE:public TET_SIM_EXAMPLE<T>
{
public:
    bool use_plasticity,use_control;
    DIAGONALIZED_PLASTICITY_MODEL_3D<T>* plasticity_model;
    T yield_ratio,plastic_clamp_ratio;
    T side_length;int resolution;
    GRID<TV> cube_grid;
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    ARRAY<VECTOR_3D<T> > plastic_goal;

    PLASTIC_CUBE_EXAMPLE()
        :use_plasticity(true),use_control(false),plasticity_model(0),yield_ratio(2),plastic_clamp_ratio(2),side_length(.5),resolution(20),
        cube_grid(resolution,resolution,resolution,-side_length/2,side_length/2,-side_length/2,side_length/2,-side_length/2,side_length/2),
        initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(5,2,0)
    {
        final_time=5;
        frame_rate=48;
        restart_step_number=0;   
        cfl_number=(T)10;
        cg_tolerance=(T)1e-2;
        use_masses_and_springs=false;use_altitude_springs=false;
        use_diagonalized_fvm=true;
        use_linear_elasticity=false;use_neo_hookean=false;
        youngs_modulus=40000;poissons_ratio=(T).45;Rayleigh_coefficient=(T).01;
        use_plasticity=true;
        use_control=true;
        yield_ratio=(T)1.2;
        plastic_clamp_ratio=3;
        strcpy(output_directory,"Plastic/output");
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=true;
        solids_parameters.collide_with_interior=false;
    }

    ~PLASTIC_CUBE_EXAMPLE()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    HEAVY_PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;

    tetrahedralized_volume.particles.Update_Position();
    tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(cube_grid);
    std::cout << "total vertices = " << particles.Size() << std::endl;std::cout << "total tetrhedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    tetrahedralized_volume.particles.Update_Position_And_Velocity();tetrahedralized_volume.particles.Store_Mass();

    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    
    if(use_control){
        plastic_goal.Resize(1,particles.Size());
        T cube_volume=tetrahedralized_volume.Total_Volume();
        for(int p=0;p<particles.Size();p++){
            T scale=particles.X(p).Lp_Norm(6)/(particles.X(p).Magnitude()+(T)1e-10);
            plastic_goal(p)=scale*particles.X(p);}
        ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::Exchange(particles.X,plastic_goal.array);
        T new_volume=tetrahedralized_volume.Total_Volume(),scale=pow(cube_volume/new_volume,(T)one_third);
        ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::Exchange(particles.X,plastic_goal.array);
        for(int p=0;p<particles.Size();p++)plastic_goal(p)*=scale;}

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<particles.Size();i++){
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
    if(!use_plasticity)return;
    if(use_control){plasticity_model=new DIAGONALIZED_PLASTICITY_CONTROL_3D<T>(*strain,yield_ratio,plastic_goal);plastic_goal.Clean_Memory();}
    else plasticity_model=new DIAGONALIZED_SIMPLE_PLASTICITY_3D<T>(*strain,yield_ratio,plastic_clamp_ratio);
    diagonalized_fvm->plasticity_model=plasticity_model;
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
