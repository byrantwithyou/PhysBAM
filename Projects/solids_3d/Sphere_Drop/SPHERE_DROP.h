//#####################################################################
// Copyright 2002, 2003, Robert Bridson, Ronald Fedkiw, Neil Molino, Zhaosheng Bao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE_DROP
//##################################################################### 
// Hanging curtain with a moving ball
#ifndef __SPHERE_DROP__
#define __SPHERE_DROP__

#include <Solids_And_Fluids/SOLIDS_3D_EXAMPLE.h>
namespace PhysBAM{

template<class T,class RW>
class SPHERE_DROP:public SOLIDS_3D_EXAMPLE<T,RW>
{
public:
    T height;
    
    SPHERE_DROP()
        :height(1.5)
    {
        last_frame=7*30;
        cfl=(T).3;
        cg_tolerance=(T)1e-3;
        cg_iterations=2000;
        
        output_directory="Sphere_Drop/output";
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=true;
        gravity=10;
    }
  
    ~SPHERE_DROP()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Surface(); 
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).triangulated_surface;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    
    std::fstream input;char filename[256];sprintf(filename,"../../Public_Data/Rigid_Bodies/sphere.tri");input.open(filename,std::ios::in|std::ios::binary);
    triangulated_surface.template Read<RW>(input);input.close();
    particles.Update_Velocity();particles.Store_Mass();
    
    triangulated_surface.Update_Bounding_Box();
    T mass_node=100;ARRAY<T>::copy(mass_node,particles.mass.array);
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.X(i)+=VECTOR_3D<T>(0,height,0);
        particles.V(i)=VECTOR_3D<T>(0,0,0);}
    
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
    SOLIDS_3D_EXAMPLE<T,RW>::Initialize_Bodies();
    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface);
    solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Linear_Finite_Volume(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface,1e6,.3,.02,true,.1,true,true);
    solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements(solids_parameters.deformable_body_parameters.list(1).triangulated_surface->triangle_mesh,1e6,1e2,true,.1);
} 
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies()
{
    int index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
};
}
#endif
