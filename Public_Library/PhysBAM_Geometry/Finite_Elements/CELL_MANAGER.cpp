//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Classes CELL_MANAGER
//#####################################################################
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_MANAGER<TV>::
CELL_MANAGER(const CELL_DOMAIN_INTERFACE<TV>& cdi_input):cdi(cdi_input)
{
    for(int s=0;s<2;s++){
        compressed[s].Resize(cdi.flat_size);
        compressed[s].Fill(-1);}
}
//#####################################################################
// Function Compress_Indices
//#####################################################################
template<class TV> void CELL_MANAGER<TV>::
Compress_Indices()
{
    if(cdi.periodic_bc)
        for(int s=0;s<2;s++){
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(cdi.grid);it.Valid();it.Next()){
                int i=cdi.Flatten(it.index);
                if(compressed[s](i)==-2) compressed[s](i)=dofs[s]++;}
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(cdi.grid,cdi.padding,GRID<TV>::GHOST_REGION,-1);it.Valid();it.Next()){
                int i=cdi.Flatten(it.index);
                if(compressed[s](i)==-2) compressed[s](i)=compressed[s](cdi.remap(i));}}
    else
        for(int s=0;s<2;s++)
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(cdi.grid,cdi.padding,GRID<TV>::WHOLE_REGION);it.Valid();it.Next()){
                int i=cdi.Flatten(it.index);
                if(compressed[s](i)==-2) compressed[s](i)=dofs[s]++;}
}
template class CELL_MANAGER<VECTOR<float,2> >;
template class CELL_MANAGER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CELL_MANAGER<VECTOR<double,2> >;
template class CELL_MANAGER<VECTOR<double,3> >;
#endif
