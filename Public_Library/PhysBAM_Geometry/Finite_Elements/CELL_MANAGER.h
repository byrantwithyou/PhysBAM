//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_MANAGER
//#####################################################################
#ifndef __CELL_MANAGER__
#define __CELL_MANAGER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV> class CELL_DOMAIN_INTERFACE;

template<class TV>
class CELL_MANAGER:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:

    const CELL_DOMAIN_INTERFACE<TV>& cdi;
    VECTOR<ARRAY<int>,2> compressed;  // inside and outside; [-1] - inactive; [-2] - active (temporary); [non-negative] - dof number (for active)
    VECTOR<int,2> dofs; // inside and outside; number of dofs

    CELL_MANAGER(const CELL_DOMAIN_INTERFACE<TV>& cdi_input);

    void Set_Active(int i,int s)
    {compressed[s](i)=-2;compressed[s](cdi.remap(i))=-2;}

    int Get_Index(TV_INT index,int s)
    {return compressed[s](cdi.Flatten(index));}

    void Compress_Indices();
};
}
#endif
