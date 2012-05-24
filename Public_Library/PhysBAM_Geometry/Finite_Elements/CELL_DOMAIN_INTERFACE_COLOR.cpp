//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_DOMAIN_INTERFACE_COLOR<TV>::
CELL_DOMAIN_INTERFACE_COLOR(const GRID<TV>& grid_input,int padding_input,int colors_input,bool wrap_input)
    :grid(grid_input),padding(padding_input),size(grid.counts+2*padding),flat_size(size.Product()),
    colors(colors_input),wrap(wrap_input)
{
    a(TV::m-1)=1;
    for(int i=TV::m-2;i>=0;i--) a(i)=a(i+1)*size(i+1);
    b=(TV_INT()+padding).Dot(a);

    remap=IDENTITY_ARRAY<>(flat_size);
    constraint_base_normal=constraint_base_tangent=0;
    for(int orientation=0;orientation<TV::m-1;orientation++){
        flat_base(orientation)=&flat_base_tangent;
        constraint_base(orientation)=&constraint_base_tangent;}
    flat_base(TV::m-1)=&flat_base_normal;
    constraint_base(TV::m-1)=&constraint_base_normal;

    if(wrap){
        for(int axis=0;axis<TV::m;axis++)
            for(int s=0;s<2;s++){
                int side=axis*2+s;
                int sign=s?-1:1;
                int diff=Flatten_Diff(sign*grid.counts(axis)*TV_INT::Axis_Vector(axis));
                for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,padding,GRID<TV>::GHOST_REGION,side);it.Valid();it.Next()){
                    int f=Flatten(it.index);
                    remap(f)=remap(f+diff);}}}

    cell_location.Resize(flat_size);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,-padding,GRID<TV>::WHOLE_REGION);it.Valid();it.Next())
        cell_location(Flatten(it.index))=-1;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,padding,GRID<TV>::GHOST_REGION,-1);it.Valid();it.Next())
        cell_location(Flatten(it.index))=1;
}
//#####################################################################
// Function Set_Flat_Base
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Set_Flat_Base_And_Resize(int extra_constraints_normal,int extra_constraints_tangent,const TV_INT& index)
{
    int flat_index=Flatten(index);

    int constraints_normal=constraint_base_normal+extra_constraints_normal;
    flat_base_normal.Resize(constraints_normal,true,true,flat_index);

    int constraints_tangent=constraint_base_tangent+extra_constraints_tangent;
    flat_base_tangent.Resize(constraints_tangent,true,true,flat_index);
}
//#####################################################################
// Function Update_Constraint_Count
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Update_Constraint_Count()
{
    constraint_base_normal=flat_base_normal.m;
    constraint_base_tangent=flat_base_tangent.m;
}
//#####################################################################
// Function Update_Total_Constraint_Count
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Update_Total_Constraint_Count()
{
    total_number_of_surface_constraints=0;
    for(int orientation=0;orientation<TV::m;orientation++)
        total_number_of_surface_constraints+=*constraint_base(orientation);
}
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >;
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >;
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >;
#endif
