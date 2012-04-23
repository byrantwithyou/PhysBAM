//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Classes CELL_MANAGER
//#####################################################################
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_DOMAIN_INTERFACE<TV>::
CELL_DOMAIN_INTERFACE(const GRID<TV>& grid_input,int padding_input,int coarse_factor_input,int interface_elements_input,int periodic_bc_input):
    grid(grid_input),padding(padding_input),size(grid.counts+2*padding),flat_size(size.Product()),coarse_factor(coarse_factor_input),
    coarse_range(TV_INT()+coarse_factor),interface_elements(interface_elements_input),periodic_bc(periodic_bc_input)
{
    a(TV::m-1)=1;
    for(int i=TV::m-2;i>=0;i--) a(i)=a(i+1)*size(i+1);
    b=(TV_INT()+padding).Dot(a);
    flat_base.Resize(interface_elements_input);
    Initialize_Remap();
}
//#####################################################################
// Function Set_Flat_Base
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE<TV>::
Set_Flat_Base(int start,int end,const TV_INT& index)
{
    int flat=Flatten(index);for(int i=start;i<end;i++) flat_base(i)=flat;
}
//#####################################################################
// Function Initialize_Remap
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE<TV>::
Initialize_Remap()
{
    remap.Resize(flat_size);
    if(periodic_bc){
        remap.Fill(-1);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
            int i=Flatten(it.index);
            remap(i)=i;}
        for(int axis=0;axis<TV::m;axis++)
            for(int s=0;s<2;s++){
                int side=axis*2+s;
                int sign=s?1:-1;
                for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,padding,GRID<TV>::GHOST_REGION,side);it.Valid();it.Next()){
                    int r=Flatten(it.index+sign*grid.counts(axis)*TV_INT::Axis_Vector(axis));
                    if(r>=0) remap(Flatten(it.index))=remap(r);}}
    }
    else for(int i=0;i<flat_size;i++) remap(i)=i;
}
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
                if(compressed[s](i)==-2) compressed[s](i)=compressed[s](cdi.Remap(i));}}
    else
        for(int s=0;s<2;s++)
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(cdi.grid,cdi.padding,GRID<TV>::WHOLE_REGION);it.Valid();it.Next()){
                int i=cdi.Flatten(it.index);
                if(compressed[s](i)==-2) compressed[s](i)=dofs[s]++;}
}
template class CELL_MANAGER<VECTOR<float,2> >;
template class CELL_MANAGER<VECTOR<float,3> >;
template class CELL_DOMAIN_INTERFACE<VECTOR<float,2> >;
template class CELL_DOMAIN_INTERFACE<VECTOR<float,3> >;
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class CELL_MANAGER<VECTOR<double,2> >;
template class CELL_MANAGER<VECTOR<double,3> >;
template class CELL_DOMAIN_INTERFACE<VECTOR<double,2> >;
template class CELL_DOMAIN_INTERFACE<VECTOR<double,3> >;
#endif
