//#####################################################################
// Copyright 2006-2007, Jeong-Mo Hong, Nipun Kwatra, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DETONATION_SHOCK_DYNAMICS
//#####################################################################
#ifndef __DETONATION_SHOCK_DYNAMICS__
#define __DETONATION_SHOCK_DYNAMICS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/GRID_AND_ARRAY_CONTAINER.h>
namespace PhysBAM{

template<class T_GRID>
class DETONATION_SHOCK_DYNAMICS
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<VECTOR<T,3> >::TYPE T_ARRAYS_RGB;
    typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
public:
    T_GRID& grid;
    const LEVELSET<TV>& levelset;
    GRID_AND_ARRAY_CONTAINER<T_GRID,T> Dn,Dn_dot,curvature,curvature_old;
    int order;

    // dsd variables that should be set before codes run
    T Dcj;
    T Dcj_min_clamp,Dcj_max_clamp;
    T A_coeff,B_coeff,C_coeff,D_coeff;
    T mutheta,dtheta;
    bool use_log_Lcj;//use true.
    int nb_width;// narrow band width to update DSD (DSD requires nb information only.)

    // Narrowband Indices
    ARRAY<TV_INT> indices_interface,indices_interface_ghost;

    BOUNDARY<TV,T> *boundary,boundary_default;
    BOUNDARY<TV,TV> *boundary_vector,boundary_vector_default;

    DETONATION_SHOCK_DYNAMICS(T_GRID& grid_input,const LEVELSET<TV>& levelset_input,const int order_input=3);
    virtual ~DETONATION_SHOCK_DYNAMICS();

    void Set_Custom_Boundary(BOUNDARY<TV,T>* boundary_input,BOUNDARY<TV,TV>* boundary_vector_input)
    {boundary=boundary_input;boundary_vector=boundary_vector_input;
    Dn.Set_Custom_Boundary(*boundary_input);Dn_dot.Set_Custom_Boundary(*boundary_input);curvature.Set_Custom_Boundary(*boundary_input);curvature_old.Set_Custom_Boundary(*boundary_input);}

//#####################################################################
    void Initialize_Grid();
    void Advance_One_Time_Step(const T_FACE_ARRAYS_SCALAR& V,const T dt,const T time,const int number_of_ghost_cells);
    void Make_NB_Indices(T_GRID &grid,T_ARRAYS_SCALAR &phi,ARRAY<TV_INT>& indices_interface,const T dt,const T time,int number_of_ghost_cells);
    bool Closest_Point_On_Boundary(T_ARRAYS_SCALAR &phi_ghost,ARRAY<TV,TV_INT> &normals_ghost,const TV& location,TV& new_location,const T tolerance=0,const int max_iterations=1) const;
    T Normal_Flame_Speed(const int axis,const TV_INT& face_index) const;
//#####################################################################
};
}
#endif
