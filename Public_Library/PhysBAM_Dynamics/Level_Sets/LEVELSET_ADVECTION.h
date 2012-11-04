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
#include <PhysBAM_Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>

namespace PhysBAM {

template<class T_GRID> struct INTERPOLATION_POLICY;
template<class T_GRID> struct BOUNDARY_POLICY;
template<class T_GRID,class T2,class T_NESTED_ADVECTION> class ADVECTION_MACCORMACK_UNIFORM;

template<class TV>
class LEVELSET_ADVECTION
{
public:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
private:
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename ADVECTION_COLLIDABLE_POLICY<GRID<TV> >::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP T_FACE_LOOKUP;typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<TV> > T_FACE_LOOKUP_COLLIDABLE;
    typedef typename REBIND<ARRAY<T,FACE_INDEX<TV::m> >,bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T,GRID<TV> > T_LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR;
public:
    LEVELSET<TV>* levelset;

    ADVECTION<GRID<TV>,T>* advection;

    T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL* nested_semi_lagrangian_collidable;
    ADVECTION_WRAPPER_COLLIDABLE_CELL<GRID<TV>,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>* semi_lagrangian_collidable;

    int reinitialization_runge_kutta_order;
    T reinitialization_cfl;
    int reinitialization_spatial_order;

    ADVECTION_MACCORMACK_UNIFORM<GRID<TV>,T,ADVECTION<GRID<TV>,T> >* advection_maccormack;

    int local_advection_spatial_order;
    bool local_semi_lagrangian_advection;

    LEVELSET_ADVECTION(LEVELSET<TV>* levelset=0);
    ~LEVELSET_ADVECTION();

    void Set_Custom_Advection(ADVECTION<GRID<TV>,T>& advection_input)
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

    void Use_Maccormack_Advection(const ARRAY<bool,TV_INT>& cell_mask);
    T Approximate_Negative_Material(const T interface_thickness=3,const T time=0) const;
    T Approximate_Positive_Material(const T interface_thickness=3,const T time=0) const;
    void Euler_Step(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocity,const T dt,const T time,const int number_of_ghost_cells);

    void Use_Local_WENO_For_Advection()
    {local_advection_spatial_order=5;local_semi_lagrangian_advection=false;}

    void Use_Local_ENO_For_Advection(const int order=3)
    {local_advection_spatial_order=order;local_semi_lagrangian_advection=false;assert(order >= 1 && order <= 3);}
            
    void Use_Local_Semi_Lagrangian_Advection()
    {local_semi_lagrangian_advection=true;local_advection_spatial_order=0;}

    void Use_Level_Set_Advection_Method()
    {local_semi_lagrangian_advection=false;local_advection_spatial_order=0;}

    void Euler_Step(const ARRAY<TV,TV_INT>& velocity,const T dt,const T time,const int number_of_ghost_cells);
    void Reinitialize(const int time_steps=10,const T time=0)
    {PhysBAM::Reinitialize(*levelset,time_steps,time,levelset->half_band_width,levelset->grid.dX.Max()*(1+min(3,local_advection_spatial_order)),reinitialization_cfl,reinitialization_runge_kutta_order,reinitialization_spatial_order,0);}

};

}
#endif
