//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_LINEAR_EXTRAPOLATION
//#####################################################################
#ifndef __BOUNDARY_LINEAR_EXTRAPOLATION__
#define __BOUNDARY_LINEAR_EXTRAPOLATION__

#include <Core/Arrays/ARRAY.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_LINEAR_EXTRAPOLATION:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef BOUNDARY<TV,T2> BASE;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
public:
    using BASE::Find_Ghost_Regions;using BASE::Boundary;

    BOUNDARY_LINEAR_EXTRAPOLATION() = default;
    virtual ~BOUNDARY_LINEAR_EXTRAPOLATION() = default;
//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const override {} // do nothing
//#####################################################################
};
}
#endif
