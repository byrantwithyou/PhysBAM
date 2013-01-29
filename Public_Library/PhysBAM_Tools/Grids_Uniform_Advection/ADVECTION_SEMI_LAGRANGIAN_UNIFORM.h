//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_UNIFORM
//#####################################################################
#ifndef __ADVECTION_SEMI_LAGRANGIAN_UNIFORM__
#define __ADVECTION_SEMI_LAGRANGIAN_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> //  T_AVERAGING=AVERAGING_UNIFORM<T_GRID>, T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2>
class ADVECTION_SEMI_LAGRANGIAN_UNIFORM:public ADVECTION<T_GRID,T2,typename T_AVERAGING::FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
public:
    using ADVECTION<T_GRID,T2,typename T_AVERAGING::FACE_LOOKUP>::Update_Advection_Equation_Cell;
    template<class T3> struct REBIND{typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T3,T_AVERAGING,typename T_INTERPOLATION::template REBIND<T3>::TYPE> TYPE;};
    template<class T_INTERPOLATION_2> struct REBIND_INTERPOLATION{typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION_2> TYPE;};

//#####################################################################
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM();
    void Update_Advection_Equation_Node(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0);
    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0);
    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost=0,const T_FACE_LOOKUP* Z_max_ghost=0,T_FACE_ARRAYS_SCALAR* Z_min=0,T_FACE_ARRAYS_SCALAR* Z_max=0);
//#####################################################################
};
}
#endif
