//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_MANAGER_COLOR
//#####################################################################
#ifndef __CELL_MANAGER_COLOR__
#define __CELL_MANAGER_COLOR__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>

namespace PhysBAM{

template<class TV> class CELL_DOMAIN_INTERFACE_COLOR;

template<class TV>
class CELL_MANAGER_COLOR
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:

    const CELL_DOMAIN_INTERFACE_COLOR<TV>& cdi;
    ARRAY<ARRAY<int> > compressed;  // [-1] - inactive; [-2] - active (temporary); [non-negative] - dof number (for active)
    ARRAY<int> dofs; // number of dofs for every color
    ARRAY<VECTOR<int,2> > uncompressed;

    CELL_MANAGER_COLOR(const CELL_DOMAIN_INTERFACE_COLOR<TV>& cdi_input);
    CELL_MANAGER_COLOR(const CELL_MANAGER_COLOR&) = delete;
    void operator=(const CELL_MANAGER_COLOR&) = delete;

    void Set_Active(int flat_index,int color)
    {compressed(color)(flat_index)=-2;compressed(color)(cdi.remap(flat_index))=-2;}

    int Get_Index(TV_INT index,int color)
    {return compressed(color)(cdi.Flatten(index));}

    int Get_Index(int flat_index,int color)
    {return compressed(color)(flat_index);}

    void Compress_Indices();
};
}
#endif
