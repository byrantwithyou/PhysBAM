#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
using namespace PhysBAM;
template<class T,class T_GRID> typename LEVELSET<T,T_GRID>::T_LINEAR_INTERPOLATION_SCALAR LEVELSET<T,T_GRID>::interpolation_default;
template<class T,class T_GRID> typename LEVELSET<T,T_GRID>::T_LINEAR_INTERPOLATION_VECTOR LEVELSET<T,T_GRID>::normal_interpolation_default;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T_GRID> LEVELSET<T,T_GRID>::
LEVELSET()
    :levelset_callbacks(0),collision_body_list(0),face_velocities_valid_mask_current(0),clamp_phi_with_collision_bodies(true),boundary_default(*new BOUNDARY_UNIFORM<T_GRID,T>),
    collision_aware_interpolation_plus(0),collision_aware_interpolation_minus(0),collision_unaware_interpolation(0),collidable_phi_replacement_value((T)1e-5)
{
    Set_Small_Number();
    Set_Max_Time_Step();
    curvature_motion=false; // default is no curvature motion
    Initialize_FMM_Initialization_Iterative_Solver();
    boundary=&boundary_default;
    interpolation=&interpolation_default;
    curvature_interpolation=&interpolation_default;
    normal_interpolation=&normal_interpolation_default;
    secondary_interpolation=0;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T_GRID> LEVELSET<T,T_GRID>::
~LEVELSET()
{
    assert(!collision_unaware_interpolation);
    delete collision_aware_interpolation_plus;
    delete collision_aware_interpolation_minus;
    delete &boundary_default;
}
//#####################################################################
// Function Set_Collision_Body_List
//#####################################################################
template<class T,class T_GRID> void LEVELSET<T,T_GRID>::
Set_Collision_Body_List(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_body_list_input,const bool set_secondary_interpolation)
{
    collision_body_list=&collision_body_list_input;
    delete collision_aware_interpolation_plus;delete collision_aware_interpolation_minus;
    collision_aware_interpolation_plus=new T_LINEAR_INTERPOLATION_SCALAR;
    collision_aware_interpolation_minus=new LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T>(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value);
    if(set_secondary_interpolation) secondary_interpolation=collision_aware_interpolation_minus;
}
template class LEVELSET<float,GRID<VECTOR<float,1> > >;
template class LEVELSET<float,GRID<VECTOR<float,2> > >;
template class LEVELSET<float,GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET<double,GRID<VECTOR<double,1> > >;
template class LEVELSET<double,GRID<VECTOR<double,2> > >;
template class LEVELSET<double,GRID<VECTOR<double,3> > >;
#endif
