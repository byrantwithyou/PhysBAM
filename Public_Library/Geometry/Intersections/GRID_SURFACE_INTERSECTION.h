//#####################################################################
// Copyright 2020, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GRID_SURFACE_INTERSECTION__
#define __GRID_SURFACE_INTERSECTION__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Grid_Tools/Grids/EDGE_INDEX.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{


template<class TV>
struct GRID_SURFACE_INTERSECTION_DATA
{
    typedef typename TV::SCALAR T;
    // whether endpoint nodes are inside object
    VECTOR<bool,2> in;

    // element index, edge theta, barycentric coordinates
    ARRAY<TRIPLE<int,T,TV> > cut_elements;
};

// Node grid, return cut edges.
template<class T,class TV>
void Grid_Surface_Intersection(
    HASHTABLE<EDGE_INDEX<2>,GRID_SURFACE_INTERSECTION_DATA<TV> >& hash,
    const GRID<TV>& grid,const SEGMENTED_CURVE_2D<T>& surface,bool compute_inside);

template<class T,class TV>
void Grid_Surface_Intersection(
    HASHTABLE<EDGE_INDEX<3>,GRID_SURFACE_INTERSECTION_DATA<TV> >& hash,
    const GRID<TV>& grid,const TRIANGULATED_SURFACE<T>& surface,bool compute_inside);
}
#endif
