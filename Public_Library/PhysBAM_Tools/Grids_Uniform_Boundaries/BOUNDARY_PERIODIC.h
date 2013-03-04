//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_PERIODIC
//#####################################################################
#ifndef __BOUNDARY_PERIODIC__
#define __BOUNDARY_PERIODIC__

#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_PERIODIC:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef typename GRID<TV>::VECTOR_INT TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef UNIFORM_GRID_ITERATOR_NODE<TV> NODE_ITERATOR;
    typedef BOUNDARY<TV,T2> BASE;
public:
    using BASE::Find_Ghost_Regions;

    BOUNDARY_PERIODIC()
    {}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
}
#endif
