//#####################################################################
// Copyright 2005, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_ITERATOR_BASE
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Grid_Tools/Grids/GRID_ITERATOR_BASE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_ITERATOR_BASE<TV>::
GRID_ITERATOR_BASE(const GRID<TV>& grid_input)
    :grid(grid_input),number_of_regions(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_ITERATOR_BASE<TV>::
GRID_ITERATOR_BASE(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input)
    :grid(grid_input),number_of_regions(0)
{
    Add_Region(region_input);Reset();
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<class TV> void GRID_ITERATOR_BASE<TV>::
Next_Helper()
{
    index(TV::dimension-1)=region.min_corner(TV::dimension-1);
    for(int i=TV::dimension-2;i>=0;i--){
        if(index(i)<region.max_corner(i)-1){index(i)++;return;}
        index(i)=region.min_corner(i);}
    Reset(current_region+1);
}
//#####################################################################
namespace PhysBAM{
template class GRID_ITERATOR_BASE<VECTOR<float,0> >;
template class GRID_ITERATOR_BASE<VECTOR<float,1> >;
template class GRID_ITERATOR_BASE<VECTOR<float,2> >;
template class GRID_ITERATOR_BASE<VECTOR<float,3> >;
template class GRID_ITERATOR_BASE<VECTOR<double,0> >;
template class GRID_ITERATOR_BASE<VECTOR<double,1> >;
template class GRID_ITERATOR_BASE<VECTOR<double,2> >;
template class GRID_ITERATOR_BASE<VECTOR<double,3> >;
}
