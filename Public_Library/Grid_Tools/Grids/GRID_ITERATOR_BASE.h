//#####################################################################
// Copyright 2005, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_ITERATOR_BASE
//#####################################################################
// Serves as base for node and cell iterators
//#####################################################################
#ifndef __GRID_ITERATOR_BASE__
#define __GRID_ITERATOR_BASE__
#include <Core/Utilities/PHYSBAM_ATTRIBUTE.h>
#include <Grid_Tools/Grids/GRID.h>

namespace PhysBAM{

template<class TV>
class GRID_ITERATOR_BASE
{
    typedef VECTOR<int,TV::dimension> TV_INT;
protected:
    const GRID<TV>& grid;
    int number_of_regions,current_region;
    RANGE<TV_INT> region,regions[1<<TV::dimension]; // at most 2^d distinct regions for iteration
    TV_INT index;
    bool valid;

    GRID_ITERATOR_BASE(const GRID<TV>& grid_input);
    GRID_ITERATOR_BASE(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input);

    void Reset_Regions()
    {number_of_regions=0;}

    void Add_Region(const RANGE<TV_INT>& region_input)
    {assert(number_of_regions<(1<<TV::dimension));
    if(TV::dimension==0) return;
    if(region_input.Edge_Lengths().Min()>=0) regions[number_of_regions++]=region_input;}

    void Reset(const int region_index=0)
    {if(TV::dimension==0){valid=true;index=TV_INT();return;}
    valid=(region_index<number_of_regions);if(valid){current_region=region_index;region=regions[current_region];index=region.Minimum_Corner();}}

public:
    bool Valid() const
    {return valid;}

    void Next() PHYSBAM_ALWAYS_INLINE
    {if(TV::dimension==0){valid=false;return;}
    if(index(TV::dimension-1)<region.max_corner(TV::dimension-1)-1) index(TV::dimension-1)++;else Next_Helper();}

    int Flat_Index()
    {int result=index(0);for(int i=1;i<TV::dimension;i++){result=result*grid.counts(i)+index(i);}return result;}

protected:
    void Next_Helper();
};
//#####################################################################
}
#endif
