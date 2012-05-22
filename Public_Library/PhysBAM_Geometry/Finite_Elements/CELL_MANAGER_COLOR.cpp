//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Classes CELL_MANAGER_COLOR
//#####################################################################
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_MANAGER_COLOR<TV>::
CELL_MANAGER_COLOR(const CELL_DOMAIN_INTERFACE_COLOR<TV>& cdi_input):cdi(cdi_input)
{
    compressed.Resize(cdi.colors);
    dofs.Resize(cdi.colors);
    for(int c=0;c<cdi.colors;c++){
        compressed(c).Resize(cdi.flat_size);
        compressed(c).Fill(-1);}
}
//#####################################################################
// Function Compress_Indices
//#####################################################################
template<class TV> void CELL_MANAGER_COLOR<TV>::
Compress_Indices()
{
    if(cdi.wrap)
        for(int c=0;c<cdi.colors;c++){
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(cdi.grid);it.Valid();it.Next()){
                int i=cdi.Flatten(it.index);
                if(compressed(c)(i)==-2) compressed(c)(i)=dofs(c)++;}
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(cdi.grid,cdi.padding,GRID<TV>::GHOST_REGION,-1);it.Valid();it.Next()){
                int i=cdi.Flatten(it.index);
                if(compressed(c)(i)==-2) compressed(c)(i)=compressed(c)(cdi.remap(i));}}
    else
        for(int c=0;c<cdi.colors;c++)
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(cdi.grid,cdi.padding,GRID<TV>::WHOLE_REGION);it.Valid();it.Next()){
                int i=cdi.Flatten(it.index);
                if(compressed(c)(i)==-2) compressed(c)(i)=dofs(c)++;}
}
template class CELL_MANAGER_COLOR<VECTOR<float,2> >;
template class CELL_MANAGER_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CELL_MANAGER_COLOR<VECTOR<double,2> >;
template class CELL_MANAGER_COLOR<VECTOR<double,3> >;
#endif
