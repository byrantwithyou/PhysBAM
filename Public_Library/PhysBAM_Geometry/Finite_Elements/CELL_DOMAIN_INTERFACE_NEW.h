//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_DOMAIN_INTERFACE_NEW
//#####################################################################
#ifndef __CELL_DOMAIN_INTERFACE_NEW__
#define __CELL_DOMAIN_INTERFACE_NEW__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class CELL_DOMAIN_INTERFACE_NEW:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    TV_INT a;
    int b;

    ARRAY<int> cell_location; // [1] - outside, [0] - boundary (inside, within padding from boundary), [-1] - inside (strictly)

public:

    const GRID<TV>& grid;
    const int padding;
    const TV_INT size;
    const int flat_size;
    const int interface_dofs;
    const bool periodic_bc;

    ARRAY<int> remap; // maps ghost cells inside for periodic bc, identity for non-periodic bc
    ARRAY<int> flat_base; // maps cut cell index to its flattened grid index
    
    CELL_DOMAIN_INTERFACE_NEW(const GRID<TV>& grid_input,int padding_input,int interface_dofs_input,bool periodic_bc_input);

    int Flatten(const TV_INT& index) const
    {return index.Dot(a)+b;}

    int Flatten_Diff(const TV_INT& index) const
    {return index.Dot(a);}

    bool Is_Outside_Cell(int i) const
    {return cell_location(i)==1;}

    bool Is_Boundary_Cell(int i) const
    {return cell_location(i)==0;}

    bool Is_Boundary_Cut_Cell(int i) const
    {return Is_Boundary_Cell(flat_base(i));}

    void Set_Flat_Base(int i,const TV_INT& index);
};
}
#endif
