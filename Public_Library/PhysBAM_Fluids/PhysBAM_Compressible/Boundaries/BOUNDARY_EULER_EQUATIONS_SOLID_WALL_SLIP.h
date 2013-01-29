//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP
//#####################################################################
#ifndef __BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP__
#define __BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP__

#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class T_GRID> class EULER_UNIFORM;

template<class TV>
class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP:public BOUNDARY<TV,VECTOR<typename TV::SCALAR,TV::m+2> >
{
    typedef typename TV::SCALAR T;typedef typename GRID<TV>::VECTOR_INT TV_INT;typedef VECTOR<T,GRID<TV>::dimension+2> TV_DIMENSION;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,GRID<TV>::dimension> TV_SIDES;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;typedef typename T_ARRAYS_BASE::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_BASE;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef VECTOR<T,2*GRID<TV>::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*GRID<TV>::dimension> TV_FACE_VECTOR;
    typedef VECTOR<bool,2*GRID<TV>::dimension> T_FACE_VECTOR_BOOL;
    typedef VECTOR<TV_DIMENSION,2*GRID<TV>::dimension> TV_DIMENSION_FACE_VECTOR;
    enum {d=TV_DIMENSION::m};
public:
    typedef BOUNDARY<TV,TV_DIMENSION> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;using BASE::Fill_Single_Ghost_Region;using BASE::Find_Ghost_Regions;using BASE::Boundary;

    EULER_UNIFORM<GRID<TV>>* euler;
    bool attenuate_using_riemann_invariants;
    T inflow_attenuation;
    bool always_attenuate;
    T_FACE_VECTOR linear_attenuations;
    T_FACE_VECTOR_BOOL linear_attenuation_faces;
    T_FACE_VECTOR S_far_field,iL_far_field,iR_far_field;
    TV_DIMENSION_FACE_VECTOR U_far_field;

    BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP(EULER_UNIFORM<GRID<TV>>* euler_input,const T_FACE_VECTOR rho_far_field,
        const T_FACE_VECTOR p_far_field,const TV_FACE_VECTOR velocity_far_field,const T inflow_attenuation_input=(T)1,
        const TV_SIDES& constant_extrapolation=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)),
        const bool always_attenuate_input=false,const T_FACE_VECTOR linear_attenuations_input=T_FACE_VECTOR(),
        const T_FACE_VECTOR_BOOL linear_attenuation_faces_input=T_FACE_VECTOR_BOOL());

//#####################################################################
    void Attenuate_To_Far_Field_Values_Using_Riemann_Invariants(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const;
    void Attenuate_To_Far_Field_Values_Using_Characteristics(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const;
    void Attenuate_To_Far_Field_Values(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const;
    void Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_BASE& u,T_ARRAYS_DIMENSION_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Fill_Single_Ghost_Region(const GRID<TV>& grid,T_ARRAYS_DIMENSION_BASE& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition_Single_Side(const GRID<TV>& grid,T_ARRAYS_DIMENSION_BASE& u,const int side,const T time) const;
    void Apply_Boundary_Condition(const GRID<TV>& grid,T_ARRAYS_DIMENSION_BASE& u,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
