//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM  
//##################################################################### 
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM__

#include <Tools/Advection/ADVECTION.h>
#include <Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<TV>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM:public ADVECTION<TV,typename TV::SCALAR,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM<TV,T,T_FACE_LOOKUP> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_UNIFORM<TV,T,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_UNIFORM<TV,typename T_FACE_LOOKUP::NESTED_LOOKUP> averaging;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list_input);
    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM();

//#####################################################################
    void Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,ARRAY<T,FACE_INDEX<TV::m> >* Z_min,ARRAY<T,FACE_INDEX<TV::m> >* Z_max);
    void Average_To_Invalidated_Face(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_values,ARRAY<bool,FACE_INDEX<TV::m> >* faces_not_to_revalidate=0);
//#####################################################################
};
}
#endif
