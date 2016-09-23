//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Nick Rasmussen, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE
//#####################################################################
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM_DEFINITION.h>
#include <Incompressible/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_FACE.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE.h>
#include <Incompressible/Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM.h>
#include <Dynamics/Advection_Equations/ADVECTION_WRAPPER_FIRE_MULTIPHASE_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <Dynamics/Interpolation/FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INCOMPRESSIBLE<TV>::
INCOMPRESSIBLE()
    :use_force_x(false),use_force_y(false),use_force_z(false),use_force(false),use_variable_surface_tension(false),use_variable_viscosity(false),
    use_variable_vorticity_confinement(false),nonzero_viscosity(false),nonzero_surface_tension(false),
    nested_semi_lagrangian_collidable(0),semi_lagrangian_collidable(0),nested_semi_lagrangian_collidable_slip(0),semi_lagrangian_collidable_slip(0),
    nested_semi_lagrangian_fire_multiphase(0),semi_lagrangian_fire_multiphase(0),nested_nested_semi_lagrangian_fire_multiphase_collidable(0),
    nested_semi_lagrangian_fire_multiphase_collidable(0),semi_lagrangian_fire_multiphase_collidable(0)
{
    advection=0; // TODO: add a default advection, possibly semi-lagrangian
    Set_Max_Time_Step();
    gravity=TV();
    Set_Surface_Tension(0);
    Set_Viscosity(0);
    Set_Vorticity_Confinement(0);
    Set_Maximum_Implicit_Viscosity_Iterations();
    Use_Explicit_Part_Of_Implicit_Viscosity();

    // for multiphase
    Use_GFM();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INCOMPRESSIBLE<TV>::
~INCOMPRESSIBLE()
{
    delete nested_semi_lagrangian_collidable;delete semi_lagrangian_collidable;
    delete nested_semi_lagrangian_collidable_slip;delete semi_lagrangian_collidable_slip;
    delete nested_semi_lagrangian_fire_multiphase;delete semi_lagrangian_fire_multiphase;
    delete nested_nested_semi_lagrangian_fire_multiphase_collidable;delete nested_semi_lagrangian_fire_multiphase_collidable;
    delete semi_lagrangian_fire_multiphase_collidable;
}
//#####################################################################
// Function Use_Semi_Lagrangian_Collidable_Advection
//#####################################################################
template<class TV> void INCOMPRESSIBLE<TV>::
Use_Semi_Lagrangian_Collidable_Advection(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list)
{
    nested_semi_lagrangian_collidable=new T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE(body_list,valid_mask);
    semi_lagrangian_collidable=new ADVECTION_WRAPPER_COLLIDABLE_FACE<TV,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE,T_FACE_LOOKUP_COLLIDABLE>(*nested_semi_lagrangian_collidable,body_list,valid_mask);
    Set_Custom_Advection(*semi_lagrangian_collidable);
}
//#####################################################################
// Function Use_Semi_Lagrangian_Collidable_Slip_Advection
//#####################################################################
template<class TV> void INCOMPRESSIBLE<TV>::
Use_Semi_Lagrangian_Collidable_Slip_Advection(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list)
{
    nested_semi_lagrangian_collidable_slip=new T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_FACE(body_list);
    semi_lagrangian_collidable_slip=new ADVECTION_WRAPPER_COLLIDABLE_FACE<TV,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_FACE,T_FACE_LOOKUP_COLLIDABLE_SLIP>(*nested_semi_lagrangian_collidable_slip,body_list,valid_mask);
    Set_Custom_Advection(*semi_lagrangian_collidable_slip);
}
//#####################################################################
// Function Use_Semi_Lagrangian_Fire_Multiphase_Advection
//#####################################################################
template<class TV> void INCOMPRESSIBLE<TV>::
Use_Semi_Lagrangian_Fire_Multiphase_Advection(const PROJECTION_DYNAMICS_UNIFORM<TV>& projection,const LEVELSET_MULTIPLE<TV>& levelset_multiple_n_plus_one)
{
    nested_semi_lagrangian_fire_multiphase=new T_ADVECTION_SEMI_LAGRANGIAN_FIRE_MULTIPHASE();
    semi_lagrangian_fire_multiphase=new T_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE(*nested_semi_lagrangian_fire_multiphase,projection,levelset_multiple_n_plus_one);
    Set_Custom_Advection(*semi_lagrangian_fire_multiphase);
}
//#####################################################################
// Function Use_Semi_Lagrangian_Fire_Multiphase_Collidable_Advection
//#####################################################################
template<class TV> void INCOMPRESSIBLE<TV>::
Use_Semi_Lagrangian_Fire_Multiphase_Collidable_Advection(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list,const PROJECTION_DYNAMICS_UNIFORM<TV>& projection,const LEVELSET_MULTIPLE<TV>& levelset_multiple_n_plus_one)
{
    nested_nested_semi_lagrangian_fire_multiphase_collidable=new T_ADVECTION_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE(body_list,valid_mask);
    nested_semi_lagrangian_fire_multiphase_collidable=new T_NESTED_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE(*nested_nested_semi_lagrangian_fire_multiphase_collidable,
        body_list,valid_mask);
    semi_lagrangian_fire_multiphase_collidable=new T_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE(*nested_semi_lagrangian_fire_multiphase_collidable,projection,
        levelset_multiple_n_plus_one);
    Set_Custom_Advection(*semi_lagrangian_fire_multiphase_collidable);
}
//#####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(T,TV,d) \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<TV> >,LINEAR_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<TV> > >;
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(float,P(VECTOR<float,1>),2)
INSTANTIATION_HELPER(float,P(VECTOR<float,2>),2);
INSTANTIATION_HELPER(float,P(VECTOR<float,3>),3);
template class INCOMPRESSIBLE<VECTOR<float,1> >;
template class INCOMPRESSIBLE<VECTOR<float,2> >;
template class INCOMPRESSIBLE<VECTOR<float,3> >;
INSTANTIATION_HELPER(double,P(VECTOR<double,1>),2);
INSTANTIATION_HELPER(double,P(VECTOR<double,2>),2);
INSTANTIATION_HELPER(double,P(VECTOR<double,3>),3);
template class INCOMPRESSIBLE<VECTOR<double,1> >;
template class INCOMPRESSIBLE<VECTOR<double,2> >;
template class INCOMPRESSIBLE<VECTOR<double,3> >;
}
