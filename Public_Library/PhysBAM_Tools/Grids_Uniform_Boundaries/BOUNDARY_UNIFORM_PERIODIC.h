//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_UNIFORM_PERIODIC
//#####################################################################
#ifndef __BOUNDARY_UNIFORM_PERIODIC__
#define __BOUNDARY_UNIFORM_PERIODIC__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class BOUNDARY_UNIFORM_PERIODIC:public BOUNDARY_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef BOUNDARY_UNIFORM<T_GRID,T2> BASE;
public:
    using BASE::Find_Ghost_Regions;

    BOUNDARY_UNIFORM_PERIODIC()
    {}

//#####################################################################
    void Fill_Ghost_Cells(const T_GRID& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const T_GRID& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
}
#endif
