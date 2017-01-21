//#####################################################################
// Copyright 2017 Lin Huang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_ITERATOR_THREADED
//#####################################################################
#include <Grid_Tools/Grids/CELL_ITERATOR_THREADED.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_ITERATOR_THREADED<TV>::
CELL_ITERATOR_THREADED(const GRID<TV>& grid_input,const int number_of_ghost_cells)
    :GRID_ITERATOR_BASE<TV>(grid_input)
{
#ifdef USE_OPENMP
    int threads=omp_get_num_threads();
    int tid=omp_get_thread_num();
#else
    int threads=1;
    int tid=0;
#endif

    Reset_Regions();
    RANGE<TV_INT> domain(grid.Domain_Indices(number_of_ghost_cells));
    int x_cells=domain.max_corner(0)-domain.min_corner(0);
    int first=x_cells*tid/threads;
    int last=x_cells*(tid+1)/threads;
    domain.max_corner(0)=domain.min_corner(0)+last;
    domain.min_corner(0)+=first;
#pragma omp critical
    if(first<last)
        Add_Region(domain);

    Reset();
}
//#####################################################################
namespace PhysBAM{
template class CELL_ITERATOR_THREADED<VECTOR<float,0> >;
template class CELL_ITERATOR_THREADED<VECTOR<float,1> >;
template class CELL_ITERATOR_THREADED<VECTOR<float,2> >;
template class CELL_ITERATOR_THREADED<VECTOR<float,3> >;
template class CELL_ITERATOR_THREADED<VECTOR<double,0> >;
template class CELL_ITERATOR_THREADED<VECTOR<double,1> >;
template class CELL_ITERATOR_THREADED<VECTOR<double,2> >;
template class CELL_ITERATOR_THREADED<VECTOR<double,3> >;
}
