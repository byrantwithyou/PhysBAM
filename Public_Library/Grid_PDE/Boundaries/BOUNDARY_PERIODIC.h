//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_PERIODIC
//#####################################################################
#ifndef __BOUNDARY_PERIODIC__
#define __BOUNDARY_PERIODIC__

#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_PERIODIC:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef BOUNDARY<TV,T2> BASE;
public:
    using BASE::Find_Ghost_Regions;
    VECTOR<bool,TV::m> is_periodic;

    BOUNDARY_PERIODIC()
    {
        is_periodic.Fill(true);
    }

    ~BOUNDARY_PERIODIC() = default;

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,
        ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const override;
//#####################################################################
};
//#####################################################################
}
#endif
