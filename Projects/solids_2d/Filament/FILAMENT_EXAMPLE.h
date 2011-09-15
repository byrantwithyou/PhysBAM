//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FILAMENT_EXAMPLE
//#####################################################################
#ifndef __FILAMENT_EXAMPLE__
#define __FILAMENT_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Forces_And_Torques/BODY_FORCES_3D.h>
namespace PhysBAM{

template<class T,class RW>
class FILAMENT_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::solids_parameters;using BASE::output_directory;using BASE::fluids_parameters;using BASE::verbose_dt;

    int hair_particles;
    T initial_height;
    T initial_orientation;
    VECTOR_2D<T> initial_velocity;
    T initial_angular_velocity;
    solid_body_collection.Add_Force(Create_Body_Forces<T>(*deformable_object.segmented_curve));
    solid_body_collection.Add_Force(Create_Edge_Springs<T>(deformable_object.segmented_curve->segment_mesh,3000,(T).9));
    //solid_body_collection.Add_Force(Create_Altitude_Springs(deformable_object.triangulated_area->triangle_mesh,300));
    //solid_body_collection.Add_Force(Create_Neo_Hookean_Elasticity(*deformable_object.triangulated_area,(T)5e4));
    //solid_body_collection.Add_Force(Create_Linear_Finite_Volume(*deformable_object.triangulated_area,3e4));
}
//#####################################################################
// Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_2D<T> >& V,const T time)
{
    V(1)=VECTOR_2D<T>(0,0);
}
//#####################################################################
// Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_2D<T> >& V,const T time)
{
    V(1)=VECTOR_2D<T>(0,0);
}
//#####################################################################
};
}
#endif
