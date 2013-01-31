//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MAC_GRID_PERIODIC
//#####################################################################
#ifndef __BOUNDARY_MAC_GRID_PERIODIC__
#define __BOUNDARY_MAC_GRID_PERIODIC__

#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_MAC_GRID_PERIODIC:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef typename GRID<TV>::VECTOR_INT TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef ARRAY<T2,FACE_INDEX<TV::m> > T_FACE_ARRAYS_T2;
    typedef UNIFORM_GRID_ITERATOR_NODE<TV> NODE_ITERATOR;
public:
    using BOUNDARY<TV,T2>::Find_Ghost_Regions;

    BOUNDARY_MAC_GRID_PERIODIC() 
    {}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Fill_Ghost_Faces(const GRID<TV>& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) PHYSBAM_OVERRIDE {} // do nothing
//#####################################################################
};
}
#endif
