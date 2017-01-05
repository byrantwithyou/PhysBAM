//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_ITERATOR_THREADED
//##################################################################### 
#include <Core/Log/LOG.h>
#include <Grid_Tools/Grids/FACE_ITERATOR_THREADED.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
// axis_input==0 means iterate through faces in all dimensions
template<class TV> FACE_ITERATOR_THREADED<TV>::
FACE_ITERATOR_THREADED(const GRID<TV>& grid_input,int number_of_ghost_cells,T_REGION region_type)
    :GRID_ITERATOR_BASE<TV>(grid_input)
{
    PHYSBAM_ASSERT(region_type==GRID<TV>::WHOLE_REGION || region_type==GRID<TV>::INTERIOR_REGION);
    int threads=omp_get_num_threads();
    int tid=omp_get_thread_num();
    
    Reset_Regions();
    RANGE<TV_INT> domain(grid.Cell_Indices(number_of_ghost_cells));
    int x_cells=domain.max_corner.x-domain.min_corner.x;

    int first=TV::m*x_cells*tid/threads;
    int last=TV::m*x_cells*(tid+1)/threads;
    for(int a=0;a<TV::m;a++){
        if(first<x_cells){
            RANGE<TV_INT> axis_domain(domain);
            axis_domain.min_corner(0)=domain.min_corner(0)+first;
            if(region_type==GRID<TV>::WHOLE_REGION) axis_domain.max_corner(a)++;
            else if(!first) axis_domain.min_corner(a)++;
            if(last<x_cells) axis_domain.max_corner(0)=domain.min_corner(0)+last;
            axes[number_of_regions]=a;
#pragma omp critical
            if(axis_domain.min_corner(0)<axis_domain.max_corner(0))
                Add_Region(axis_domain);}
        last-=x_cells;
        if(last<0) break;
        first-=x_cells;
        if(first<0) first=0;}
    axis=axes[0];
    Reset();
}
//#####################################################################
namespace PhysBAM{
template class FACE_ITERATOR_THREADED<VECTOR<float,1> >;
template class FACE_ITERATOR_THREADED<VECTOR<float,2> >;
template class FACE_ITERATOR_THREADED<VECTOR<float,3> >;
template class FACE_ITERATOR_THREADED<VECTOR<double,1> >;
template class FACE_ITERATOR_THREADED<VECTOR<double,2> >;
template class FACE_ITERATOR_THREADED<VECTOR<double,3> >;
}
