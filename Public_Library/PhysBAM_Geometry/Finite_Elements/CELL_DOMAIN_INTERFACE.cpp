//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_DOMAIN_INTERFACE<TV>::
CELL_DOMAIN_INTERFACE(const GRID<TV>& grid_input,int padding_input,int coarse_factor_input,int interface_elements_input,bool periodic_bc_input)
    :grid(grid_input),padding(padding_input),size(grid.counts+2*padding),flat_size(size.Product()),coarse_factor(coarse_factor_input),
    interface_elements(interface_elements_input),periodic_bc(periodic_bc_input),coarse_range(TV_INT(),TV_INT()+coarse_factor)
{
    a(TV::m-1)=1;
    for(int i=TV::m-2;i>=0;i--) a(i)=a(i+1)*size(i+1);
    b=(TV_INT()+padding).Dot(a);
    flat_coarse_offset=Flatten_Diff(TV_INT()+coarse_factor-1);
    flat_base.Resize(interface_elements_input);
    Initialize();
}
//#####################################################################
// Function Set_Flat_Base
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE<TV>::
Set_Flat_Base(int start,int end,const TV_INT& index)
{
    int flat=Flatten(index);
    bool bdy=!(Is_Inside_Cell(flat) && Is_Inside_Cell(flat+flat_coarse_offset));
    for(int i=start;i<end;i++){
        flat_base(i)=flat;
        bdy_element(i)=bdy;}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE<TV>::
Initialize()
{
    remap=IDENTITY_ARRAY<>(flat_size);

    if(periodic_bc){
        for(int axis=0;axis<TV::m;axis++)
            for(int s=0;s<2;s++){
                int side=axis*2+s;
                int sign=s?-1:1;
                int s=Flatten_Diff(sign*grid.counts(axis)*TV_INT::Axis_Vector(axis));
                for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,padding,GRID<TV>::GHOST_REGION,side);it.Valid();it.Next()){
                    int f=Flatten(it.index);
                    remap(f)=remap(f+s);}}}

    cell_location.Resize(flat_size);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,-padding,GRID<TV>::WHOLE_REGION);it.Valid();it.Next())
        cell_location(Flatten(it.index))=-1;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,padding,GRID<TV>::GHOST_REGION,-1);it.Valid();it.Next())
        cell_location(Flatten(it.index))=1;
    bdy_element.Resize(interface_elements);
}
template class CELL_DOMAIN_INTERFACE<VECTOR<float,2> >;
template class CELL_DOMAIN_INTERFACE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CELL_DOMAIN_INTERFACE<VECTOR<double,2> >;
template class CELL_DOMAIN_INTERFACE<VECTOR<double,3> >;
#endif
