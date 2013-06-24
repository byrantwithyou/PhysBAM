//#####################################################################
// Copyright 2003-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Geometry/Spatial_Acceleration/BOX_PARTITION.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BOX_PARTITION<TV>::
BOX_PARTITION(const RANGE<TV>& box,const TV_INT& counts)
    :grid(counts,box,true),cells(grid.Domain_Indices())
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BOX_PARTITION<TV>::
~BOX_PARTITION()
{
}
//#####################################################################
// Function Add
//#####################################################################
template<class TV> void BOX_PARTITION<TV>::
Add(const RANGE<TV>& range,int id)
{
    for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices().Intersect(Range(range)));it.Valid();it.Next())
        cells(it.index).Set(id);
}
//#####################################################################
// Function Remove
//#####################################################################
template<class TV> void BOX_PARTITION<TV>::
Remove(const RANGE<TV>& range,int id)
{
    for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices().Intersect(Range(range)));it.Valid();it.Next())
        cells(it.index).Delete_If_Present(id);
}
//#####################################################################
// Function Find
//#####################################################################
template<class TV> void BOX_PARTITION<TV>::
Find(const RANGE<TV>& range,ARRAY<int>& ids) const
{
    for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices().Intersect(Range(range)));it.Valid();it.Next())
        cells(it.index).Append_Keys(ids);
    ids.Prune_Duplicates();
}
//#####################################################################
// Function Find
//#####################################################################
template<class TV> void BOX_PARTITION<TV>::
Find(const TV& X,ARRAY<int>& ids) const
{
    cells(grid.Clamp_To_Cell(X)).Append_Keys(ids);
}
//#####################################################################
// Function Remove_All
//#####################################################################
template<class TV> void BOX_PARTITION<TV>::
Remove_All()
{
    for(int i=0;i<cells.array.m;i++)
        cells.array(i).Clean_Memory();
}
template class BOX_PARTITION<VECTOR<float,1> >;
template class BOX_PARTITION<VECTOR<float,2> >;
template class BOX_PARTITION<VECTOR<float,3> >;
template class BOX_PARTITION<VECTOR<double,1> >;
template class BOX_PARTITION<VECTOR<double,2> >;
template class BOX_PARTITION<VECTOR<double,3> >;
}
