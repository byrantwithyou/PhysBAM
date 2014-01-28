//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_DOMAIN_INTERFACE_COLOR<TV>::
CELL_DOMAIN_INTERFACE_COLOR(const GRID<TV>& grid_input,int padding_input,int colors_input)
    :grid(grid_input),padding(padding_input),size(grid.counts+2*padding),flat_size(size.Product()),
    colors(colors_input)
{
    a(TV::m-1)=1;
    for(int i=TV::m-2;i>=0;i--) a(i)=a(i+1)*size(i+1);
    b=(TV_INT()+padding).Dot(a);

    remap=IDENTITY_ARRAY<>(flat_size);
    constraint_base_n=constraint_base_t=constraint_base_scalar=0;
    nc_present=dc_present=sc_present=false;
    for(int orientation=0;orientation<TV::m-1;orientation++){
        flat_base(orientation)=&flat_base_t;
        constraint_base(orientation)=&constraint_base_t;}
    flat_base(TV::m-1)=&flat_base_n;
    constraint_base(TV::m-1)=&constraint_base_n;

    for(int axis=0;axis<TV::m;axis++)
        for(int s=0;s<2;s++){
            int side=axis*2+s;
            int sign=s?-1:1;
            int diff=Flatten_Diff(sign*grid.counts(axis)*TV_INT::Axis_Vector(axis));
            for(CELL_ITERATOR<TV> it(grid,padding,GRID<TV>::GHOST_REGION,side);it.Valid();it.Next()){
                int f=Flatten(it.index);
                remap(f)=remap(f+diff);}}

    cell_location.Resize(flat_size);
    for(CELL_ITERATOR<TV> it(grid,-padding,GRID<TV>::WHOLE_REGION);it.Valid();it.Next())
        cell_location(Flatten(it.index))=-1;
    for(CELL_ITERATOR<TV> it(grid,padding,GRID<TV>::GHOST_REGION,-1);it.Valid();it.Next())
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
// Function Construct_Surface_Meshes
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Construct_Surface_Meshes(const GRID<TV>& phi_grid,const ARRAY<T,TV_INT>& phi_value,const ARRAY<int,TV_INT>& phi_color)
{
    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();
    MARCHING_CUBES_COLOR<TV>::Get_Elements(index_to_cell_elements,phi_grid,phi_color,phi_value,phi_grid.counts.Product());
}
//#####################################################################
// Function Interpolate_Level_Set_To_Double_Fine_Grid
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Interpolate_Level_Set_To_Double_Fine_Grid(const RANGE<TV_INT>& range_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,
    const RANGE<TV_INT>& range,ARRAY<T,TV_INT>& phi_value,ARRAY<int,TV_INT>& phi_color,T threshold)
{
    PHYSBAM_ASSERT(range_input.Edge_Lengths()*2==range.Edge_Lengths()+1);
    TV_INT copy_offset=range.min_corner-range_input.min_corner*2;
    for(RANGE_ITERATOR<TV_INT::m> it(range_input);it.Valid();it.Next()){
        phi_value(it.index*2+copy_offset)=phi_value_input(it.index);
        phi_color(it.index*2+copy_offset)=phi_color_input(it.index);}

    RANGE<TV_INT> source_range(TV_INT(),range_input.Edge_Lengths());
    TV_INT scale(TV_INT()+2);
    for(int axis=0;axis<TV_INT::m;axis++){
        source_range.max_corner(axis)--;
        for(RANGE_ITERATOR<TV_INT::m> it(source_range);it.Valid();it.Next()){
            TV_INT a=scale*it.index+range.min_corner,m(a),b(a);
            m(axis)++;
            b(axis)+=2;
            T phi_value_a=phi_value(a);
            T phi_value_b=phi_value(b);
            int phi_color_a=phi_color(a);
            int phi_color_b=phi_color(b);
            if(phi_color_a!=phi_color_b){
                phi_value(m)=(T).5*abs(phi_value_a-phi_value_b);
                phi_color(m)=(phi_value_a<phi_value_b)?phi_color_b:phi_color_a;}
            else{
                phi_value(m)=(T).5*abs(phi_value_a+phi_value_b);
                phi_color(m)=phi_color_a;}}
        source_range.max_corner(axis)=range.Edge_Lengths()(axis);
        scale(axis)=1;}

    // Perturb levelset
    for(RANGE_ITERATOR<TV_INT::m> it(phi_value.domain);it.Valid();it.Next()){
        T& value=phi_value(it.index);
        if(value<threshold) value=threshold;}
}
//#####################################################################
// Function Interpolate_Level_Set_To_Double_Fine_Grid
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Interpolate_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,
    const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,
    const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,ARRAY<int,TV_INT>& phi_color,T tol)
{
    Interpolate_Level_Set_To_Double_Fine_Grid(phi_grid_input.Node_Indices(),phi_value_input,
        phi_color_input,phi_grid.Node_Indices(),phi_value,phi_color,(phi_grid.dX.Min()/phi_grid.counts.Max())*tol);
}
//#####################################################################
// Function Interpolate_Level_Set_To_Double_Fine_Grid
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Interpolate_Mac_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,
    const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,
    const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,ARRAY<int,TV_INT>& phi_color,T tol)
{
    PHYSBAM_ASSERT(phi_value_input.domain.min_corner(0)<=-1); // Require a layer of ghost cells
    PHYSBAM_ASSERT(phi_grid_input.Is_MAC_Grid());
    Interpolate_Level_Set_To_Double_Fine_Grid(phi_grid_input.Domain_Indices(1),phi_value_input,
        phi_color_input,phi_grid.Node_Indices(1),phi_value,phi_color,(phi_grid.dX.Min()/phi_grid.counts.Max())*tol);
}
//#####################################################################
// Function Interpolate_Level_Set_To_Double_Fine_Grid
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Interpolate_Level_Set_To_Double_Fine_Grid(const RANGE<TV_INT>& range_input,
    const ARRAY<T,TV_INT>& phi_value_input,const RANGE<TV_INT>& range,ARRAY<T,TV_INT>& phi_value,T threshold)
{
    PHYSBAM_ASSERT(range_input.Edge_Lengths()*2==range.Edge_Lengths()+1);
    TV_INT copy_offset=range.min_corner-range_input.min_corner*2;
    for(RANGE_ITERATOR<TV_INT::m> it(range_input);it.Valid();it.Next())
        phi_value(it.index*2+copy_offset)=phi_value_input(it.index);

    RANGE<TV_INT> source_range(TV_INT(),range_input.Edge_Lengths());
    TV_INT scale(TV_INT()+2);
    for(int axis=0;axis<TV_INT::m;axis++){
        source_range.max_corner(axis)--;
        for(RANGE_ITERATOR<TV_INT::m> it(source_range);it.Valid();it.Next()){
            TV_INT a=scale*it.index+range.min_corner,m(a),b(a);
            m(axis)++;
            b(axis)+=2;
            phi_value(m)=(T).5*(phi_value(a)+phi_value(b));}
        source_range.max_corner(axis)=range.Edge_Lengths()(axis);
        scale(axis)=1;}

    // Perturb levelset
    for(RANGE_ITERATOR<TV_INT::m> it(phi_value.domain);it.Valid();it.Next()){
        T& value=phi_value(it.index);
        if(abs(value)<threshold) value=sign_nonzero(value)*threshold;}
}
//#####################################################################
// Function Interpolate_Level_Set_To_Double_Fine_Grid
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Interpolate_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,const ARRAY<T,TV_INT>& phi_value_input,const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,T tol)
{
    Interpolate_Level_Set_To_Double_Fine_Grid(phi_grid_input.Node_Indices(),phi_value_input,phi_grid.Node_Indices(),phi_value,(phi_grid.dX.Min()/phi_grid.counts.Max())*tol);
}
//#####################################################################
// Function Interpolate_Level_Set_To_Double_Fine_Grid
//#####################################################################
template<class TV> void CELL_DOMAIN_INTERFACE_COLOR<TV>::
Interpolate_Mac_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,const ARRAY<T,TV_INT>& phi_value_input,const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,T tol)
{
    PHYSBAM_ASSERT(phi_value_input.domain.min_corner(0)<=-1); // Require a layer of ghost cells
    PHYSBAM_ASSERT(phi_grid_input.Is_MAC_Grid());
    Interpolate_Level_Set_To_Double_Fine_Grid(phi_grid_input.Domain_Indices(1),phi_value_input,phi_grid.Node_Indices(1),phi_value,(phi_grid.dX.Min()/phi_grid.counts.Max())*tol);
}
namespace PhysBAM{
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >;
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >;
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >;
template class CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >;
}
