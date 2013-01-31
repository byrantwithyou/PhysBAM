//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_REFLECTION_UNIFORM
//#####################################################################
#ifndef __BOUNDARY_REFLECTION_UNIFORM__
#define __BOUNDARY_REFLECTION_UNIFORM__

#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_REFLECTION_UNIFORM:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef typename GRID<TV>::VECTOR_INT TV_INT;typedef VECTOR<bool,GRID<TV>::dimension> TV_BOOL;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,GRID<TV>::dimension> TV_SIDES;typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef UNIFORM_GRID_ITERATOR_NODE<TV> NODE_ITERATOR;
public:
    typedef BOUNDARY<TV,T2> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;using BASE::Fill_Single_Ghost_Region;using BASE::Find_Ghost_Regions;
    using BASE::Boundary;

    BOUNDARY_REFLECTION_UNIFORM(const TV_SIDES& constant_extrapolation=TV_SIDES())
    {
        Set_Constant_Extrapolation(constant_extrapolation);
    }

//#####################################################################
    virtual void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    virtual void Fill_Single_Ghost_Region(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells=3) const PHYSBAM_OVERRIDE ;
    virtual void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) PHYSBAM_OVERRIDE {} // do nothing
//#####################################################################
};
}
#endif
