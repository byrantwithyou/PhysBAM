//#####################################################################
// Copyright 2003-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOX_PARTITION
//#####################################################################
#ifndef __BOX_PARTITION__
#define __BOX_PARTITION__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/RANGE.h>
#include <Grid_Tools/Grids/GRID.h>
namespace PhysBAM{

template<class TV>
class BOX_PARTITION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    GRID<TV> grid;
    ARRAY<HASHTABLE<int>,TV_INT> cells;

    BOX_PARTITION(const RANGE<TV>& box,const TV_INT& counts);
    ~BOX_PARTITION();
    void Add(const RANGE<TV>& range,int id);
    void Remove(const RANGE<TV>& range,int id);
    void Find(const RANGE<TV>& range,ARRAY<int>& ids) const;
    void Find(const TV& X,ARRAY<int>& ids) const;
    void Remove_All();

    RANGE<TV_INT> Range(const RANGE<TV>& box) const
    {return grid.Clamp_To_Cell(box);}

//#####################################################################
};
}
#endif
