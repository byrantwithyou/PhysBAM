//#####################################################################
// Copyright 2003-2006, Ron Fedkiw, Neil Molino, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TORUS_EXAMPLE
//#####################################################################
#ifndef __EMBEDDED_TORUS_EXAMPLE__
#define __EMBEDDED_TORUS_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "ANALYTIC_TORUS.h"
#include <Forces_And_Torques/BODY_FORCES_3D.h>
namespace PhysBAM{

template <class T,class RW>
class EMBEDDED_TORUS_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;

    ANALYTIC_TORUS<T> torus;
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int m_input,n_input,mn_input;
    bool cube_mesh;
    bool fine_torus;

    EMBEDDED_TORUS_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
        cube_mesh(false),fine_torus(false)
    {
        m_input=13;n_input=13;mn_input=5;
        if(fine_torus){m_input=37;n_input=37;mn_input=13;}
        last_frame=240;
        restart=false;restart_frame=3;  
        output_directory="Embedded_Torus/output";

        solids_parameters.cfl=(T)10;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.check_initial_mesh_for_self_intersection=false;
        solids_parameters.use_constant_mass=true;
        solids_parameters.collision_repulsion_spring_multiplier=100;
        solids_parameters.collisions_nonrigid_collision_attempts=3;
        solids_parameters.collisions_repulsion_spring_constant_over_mass_times_length=30;
    }

    ~EMBEDDED_TORUS_EXAMPLE()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
virtual void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Tetrahedralized_Volume(15); // 15 is for the hashtable size
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(index);
    EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_tetrahedralized_volume=*deformable_object.embedded_tetrahedralized_volume;
    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& embedded_tetrahedralized_volume_boundary_surface=*deformable_object.embedded_tetrahedralized_volume_boundary_surface;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedded_tetrahedralized_volume.tetrahedralized_volume;
    PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;
    
    GRID<TV> torus_grid(m_input,n_input,mn_input,torus.bounding_box);
    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(torus_grid);else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(torus_grid);
    ARRAY<T> phi(particles.array_collection->Size());for(int p=0;p<phi.m;p++) phi(p)=torus.Phi(particles.X(p));
    embedded_tetrahedralized_volume.Calculate_Triangulated_Surface_From_Levelset_On_Tetrahedron_Nodes(phi);
       
    tetrahedralized_volume.Set_Density(500);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);   
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=1;i<=particles.array_collection->Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}
    embedded_tetrahedralized_volume_boundary_surface.Create_Boundary_Surface_From_Manifold_Embedded_Surface();

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();

    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.embedded_tetrahedralized_volume->tetrahedralized_volume;
    solid_body_collection.Add_Force(Create_Body_Forces<T>(tetrahedralized_volume));
    if(fine_torus) solid_body_collection.Add_Force(Create_Edge_Springs<T>(tetrahedralized_volume,60,2));
    else solid_body_collection.Add_Force(Create_Edge_Springs<T>(tetrahedralized_volume,600,2));
    if(fine_torus) solid_body_collection.Add_Force(Create_Altitude_Springs<T>(tetrahedralized_volume,6));
    else solid_body_collection.Add_Force(Create_Altitude_Springs<T>(tetrahedralized_volume,300));
    //deformable_object.Add_Diagonalized_Linear_Elasticity(tetrahedralized_volume);
    //deformable_object.Add_Diagonalized_Neo_Hookean_Elasticity(tetrahedralized_volume);
    //deformable_object.Add_Diagonalized_Linear_Finite_Volume(tetrahedralized_volume);    
} 
//#####################################################################
// Function Write_Data_Files
//#####################################################################
void Write_Data_Files(const T time,const int frame) const
{
    BASE::Write_Data_Files(time,frame);
    
    // Additional Outputs
    if(solids_parameters.deformable_body_parameters.list(1).embedded_tetrahedralized_volume_boundary_surface){
        FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/boundary_surface%d.tri",output_directory.c_str(),frame),
            solids_parameters.deformable_body_parameters.list(1).embedded_tetrahedralized_volume_boundary_surface->boundary_surface);}
}
//#####################################################################
// Function Read_Data_Files
//#####################################################################
void Read_Data_Files(const T time,const int frame) const
{
    BASE::Read_Data_Files(time,frame);

    // Additional Inputs
    if(solids_parameters.deformable_body_parameters.list(1).embedded_tetrahedralized_volume_boundary_surface){
        FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/boundary_surface%d.tri",output_directory.c_str(),frame),
            solids_parameters.deformable_body_parameters.list(1).embedded_tetrahedralized_volume_boundary_surface->boundary_surface);}
}
//#####################################################################
};
}
#endif
