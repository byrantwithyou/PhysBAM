//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_UNIFORM
//#####################################################################
#ifndef __ADVECTION_SEMI_LAGRANGIAN_UNIFORM__
#define __ADVECTION_SEMI_LAGRANGIAN_UNIFORM__

#include <Grid_PDE/Advection/ADVECTION.h>
#include <Grid_PDE/Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA.h>
#include <Grid_PDE/Interpolation/AVERAGING_UNIFORM.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> //  T_AVERAGING=AVERAGING_UNIFORM<TV>, T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<TV,T2>
class ADVECTION_SEMI_LAGRANGIAN_UNIFORM:public ADVECTION<TV,T2,typename T_AVERAGING::FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
public:
    using ADVECTION<TV,T2,typename T_AVERAGING::FACE_LOOKUP>::Update_Advection_Equation_Cell;
    template<class T_INTERPOLATION_2> struct REBIND_INTERPOLATION{typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T2,T_AVERAGING,T_INTERPOLATION_2> TYPE;};

//#####################################################################
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM();
    void Update_Advection_Equation_Node(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0);
    void Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0);
    void Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost=0,const T_FACE_LOOKUP* Z_max_ghost=0,ARRAY<T,FACE_INDEX<TV::m> >* Z_min=0,ARRAY<T,FACE_INDEX<TV::m> >* Z_max=0);
//#####################################################################
};
}
#endif
