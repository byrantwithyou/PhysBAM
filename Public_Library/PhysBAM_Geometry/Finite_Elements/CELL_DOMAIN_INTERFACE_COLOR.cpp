//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
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
    constraint_base_n=constraint_base_t=0;
    for(int orientation=0;orientation<TV::m-1;orientation++){
        flat_base(orientation)=&flat_base_t;
        constraint_base(orientation)=&constraint_base_t;}
    flat_base(TV::m-1)=&flat_base_n;
    constraint_base(TV::m-1)=&constraint_base_n;

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
// Function Set_Flat_Base_And_Resize
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Set_Flat_Base_And_Resize(int extra_constraints_n,int extra_constraints_t,const TV_INT& index)
{
    int flat_index=Flatten(index);
    flat_base_n.Resize(constraint_base_n+extra_constraints_n,true,true,flat_index);
    flat_base_t.Resize(constraint_base_t+extra_constraints_t,true,true,flat_index);
}
//#####################################################################
// Function Set_Flat_Base_And_Resize_Scalar
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Set_Flat_Base_And_Resize_Scalar(int extra_constraints_scalar,const TV_INT& index)
{
    flat_base_scalar.Resize(constraint_base_scalar+extra_constraints_scalar,true,true,Flatten(index));
}
//#####################################################################
// Function Update_Constraint_Count
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Update_Constraint_Count()
{
    constraint_base_n=flat_base_n.m;
    constraint_base_t=flat_base_t.m;
}
//#####################################################################
// Function Update_Constraint_Count_Scalar
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Update_Constraint_Count_Scalar()
{
    constraint_base_scalar=flat_base_scalar.m;
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
//#####################################################################
// Function Interpolate_Level_Set_To_Double_Fine_Grid
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Interpolate_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,
    const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,
    const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,ARRAY<int,TV_INT>& phi_color,T tol)
{
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(phi_grid_input);it.Valid();it.Next()){
        phi_value(it.index*2)=phi_value_input(it.index);
        phi_color(it.index*2)=phi_color_input(it.index);}
    
    TV_INT counts(phi_grid_input.counts+1);
    TV_INT scale(TV_INT()+2);
    for(int axis=0;axis<TV::m;axis++){
        counts(axis)=phi_grid_input.counts(axis);
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),counts));it.Valid();it.Next()){
            TV_INT a=scale*it.index;
            TV_INT m(a);m+=TV_INT::Axis_Vector(axis);
            TV_INT b(m);b+=TV_INT::Axis_Vector(axis);
            const T& phi_value_a=phi_value(a);
            const T& phi_value_b=phi_value(b);
            const int& phi_color_a=phi_color(a);
            const int& phi_color_b=phi_color(b);
            if(phi_color_a!=phi_color_b){
                phi_value(m)=(T).5*abs(phi_value_a-phi_value_b);
                phi_color(m)=(phi_value_a<phi_value_b)?phi_color_b:phi_color_a;}
            else{
                phi_value(m)=(T).5*abs(phi_value_a+phi_value_b);
                phi_color(m)=phi_color_a;}}
        counts(axis)=phi_grid.counts(axis)+1;
        scale(axis)=1;}
    
    // Perturb levelset
    T panic_threshold=phi_grid.dX.Min()*tol;
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(phi_grid);it.Valid();it.Next()){
        T& value=phi_value(it.index);
        if(value<panic_threshold) value=panic_threshold;}
}
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >;
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >;
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >;
#endif
