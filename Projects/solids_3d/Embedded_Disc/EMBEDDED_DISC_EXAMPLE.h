//#####################################################################
// Copyright 2003-2005, Zhaosheng Bao, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_DISC_EXAMPLE
//#####################################################################
#ifndef __EMBEDDED_DISC_EXAMPLE__
#define __EMBEDDED_DISC_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <Constitutive_Models/DIAGONALIZED_SIMPLE_PLASTICITY_2D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template <class T,class RW>
class EMBEDDED_DISC_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;

    T initial_height;
    bool use_plasticity,use_bending_plasticity;

    EMBEDDED_DISC_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE),initial_height(1)
    {
        last_frame=240;
        restart=false;restart_frame=3;  
        solids_parameters.cfl=(T).25;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Embedded_Disc/output";
        solids_parameters.check_initial_mesh_for_self_intersection=false;

        solids_parameters.use_constant_mass=true;
/*
        solids_parameters.collision_repulsion_spring_multiplier=100;
        solids_parameters.collisions_nonrigid_collision_attempts=3;      
        solids_parameters.collisions_repulsion_spring_constant_over_mass_times_length=1e-10;
*/
        solids_parameters.collisions_disable_repulsions_based_on_proximity_factor=1.5;

        solids_parameters.perform_self_collision=true;
        use_plasticity=true;
        use_bending_plasticity=true;
    }

    ~EMBEDDED_DISC_EXAMPLE()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Triangulated_Surface();
    EMBEDDED_TRIANGULATED_SURFACE<T>& embedded_triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).embedded_triangulated_surface;
    TRIANGLES_OF_MATERIAL_3D<T>& triangles_of_material=*solids_parameters.deformable_body_parameters.list(index).triangles_of_material;
    TRIANGULATED_SURFACE<T>& triangulated_surface=embedded_triangulated_surface.triangulated_surface;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    
    FILE_UTILITIES::Read_From_File<RW>(data_directory+"/Rigid_Bodies/sphere.tri",triangulated_surface);
    particles.Update_Velocity();particles.Store_Mass();
    
    triangulated_surface.Update_Bounding_Box();
    T mass_node=100;ARRAY<T>::copy(mass_node,particles.mass.array);
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.X(i)+=VECTOR_3D<T>(0,initial_height,0);
        particles.V(i)=VECTOR_3D<T>();}

    // Scoop out Disc (or cap) from Sphere
    VECTOR_3D<T>center=triangulated_surface.bounding_box->Center();
    T radius=(T).5*triangulated_surface.bounding_box->Size().Magnitude();
    center.x+=(T).5*triangulated_surface.bounding_box->Size().x;
    SPHERE<T> sphere(center,radius);
    ARRAY<T> phi(particles.array_collection->Size());
    for(int p=0;p<particles.array_collection->Size();p++)phi(p)=sphere.Signed_Distance(particles.X(p));
    embedded_triangulated_surface.Calculate_Segmented_Curve_From_Levelset_On_Triangle_Nodes(phi);
    triangles_of_material.Create_Material_Surface();
    
    // rigid bodies
    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    Get_Initial_Data();
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces((solids_parameters.deformable_body_parameters.list(1).embedded_triangulated_surface)->triangulated_surface);
    DIAGONALIZED_PLASTICITY_MODEL_2D<T>* plasticity_model=0;
    if(use_plasticity) plasticity_model=new DIAGONALIZED_SIMPLE_PLASTICITY_2D<T>((solids_parameters.deformable_body_parameters.list(1).embedded_triangulated_surface->triangulated_surface).triangle_mesh.triangles.m,2,2);
    solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Linear_Finite_Volume(solids_parameters.deformable_body_parameters.list(1).embedded_triangulated_surface->triangulated_surface,1e6,.3,.02,true,.1,true,true,plasticity_model);
    solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements((solids_parameters.deformable_body_parameters.list(1).embedded_triangulated_surface->triangulated_surface).triangle_mesh,1e5,1e3,true,.1,use_bending_plasticity,.03,.03,0,0);
} 
//#####################################################################
};
}
#endif
