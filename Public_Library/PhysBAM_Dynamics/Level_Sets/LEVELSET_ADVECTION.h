//#####################################################################
// Copyright 2009, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_ADVECTION
//##################################################################### 
#ifndef __LEVELSET_ADVECTION__
#define __LEVELSET_ADVECTION__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Geometry/Advection_Collidable/ADVECTION_COLLIDABLE_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>

namespace PhysBAM {

template<class T_GRID> struct INTERPOLATION_POLICY;
template<class T_GRID> struct BOUNDARY_POLICY;

template<class T_GRID>
class LEVELSET_ADVECTION
{
public:
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
private:
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> T_FACE_LOOKUP_COLLIDABLE;
    typedef typename REBIND<ARRAY<T,FACE_INDEX<TV::m> >,bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T,T_GRID> T_LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR;
public:
    LEVELSET<TV>* levelset;

    ADVECTION<T_GRID,T>* advection;

    T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL* nested_semi_lagrangian_collidable;
    ADVECTION_WRAPPER_COLLIDABLE_CELL<T_GRID,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>* semi_lagrangian_collidable;

    int reinitialization_runge_kutta_order;
    T reinitialization_cfl;
    int reinitialization_spatial_order;

    LEVELSET_ADVECTION(LEVELSET<TV>* _levelset);
    ~LEVELSET_ADVECTION();

    void Set_Custom_Advection(ADVECTION<T_GRID,T>& advection_input)
    {advection=&advection_input;}
    void Use_Semi_Lagrangian_Collidable_Advection(const T_GRID_BASED_COLLISION_GEOMETRY& body_list,const T phi_replacement_value,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input);
    void HJ_WENO(const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus) const;
    void HJ_ENO(const int order,const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus) const; 

    void Set_Reinitialization_Runge_Kutta_Order(const int order=3)
    {assert(order >=1 && order <=3);reinitialization_runge_kutta_order=order;}

    void Set_Reinitialization_CFL(const T cfl=.5)
    {reinitialization_cfl=cfl;assert(cfl <= 1);}

    void Use_WENO_For_Reinitialization() // 5th order
    {reinitialization_spatial_order=5;}

    void Use_ENO_For_Reinitialization(const int order=3)
    {assert(order >=1 && order <= 3);reinitialization_spatial_order=order;}
};

}
#endif
