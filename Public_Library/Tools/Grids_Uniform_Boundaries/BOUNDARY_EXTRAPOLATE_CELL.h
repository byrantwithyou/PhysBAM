//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EXTRAPOLATE_CELL
//#####################################################################
#ifndef __BOUNDARY_EXTRAPOLATE_CELL__
#define __BOUNDARY_EXTRAPOLATE_CELL__

#include <Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_EXTRAPOLATE_CELL:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:

    BOUNDARY_EXTRAPOLATE_CELL()
    {}

    virtual ~BOUNDARY_EXTRAPOLATE_CELL()
    {}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const override;
//#####################################################################
};
//#####################################################################
}
#endif
