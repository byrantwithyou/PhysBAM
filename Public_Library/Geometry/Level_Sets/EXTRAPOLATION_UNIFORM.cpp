//####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_UNIFORM  
//##################################################################### 
#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_UNIFORM<TV,T2>::
EXTRAPOLATION_UNIFORM(const GRID<TV>& grid,const T_ARRAYS_BASE& phi_input,ARRAYS_ND_BASE<T2,TV_INT>& u_input,const int ghost_cells_input)
    :u(u_input),phi(phi_input),seed_indices(0),seed_done(0),collision_aware_extrapolation(false),neighbors_visible(0),ghost_cells(ghost_cells_input)
{
    Set_Band_Width();
    Set_Isobaric_Fix_Width();
    Set_Small_Number();
    node_grid=grid.Is_MAC_Grid()?grid.Get_Regular_Grid_At_MAC_Positions():grid;
    TV DX2=node_grid.dX*node_grid.dX;
    if(TV::m>=2){optimization_scale[2]=DX2[1]/DX2[0]; // dy^2/dx^2
        if(TV::m==3){optimization_scale[0]=DX2[2]/DX2[1];optimization_scale[1]=DX2[2]/DX2[0];}} // dz^2/dy^2 and dz^2/dx^2 
    RANGE<TV_INT> domain_indices=node_grid.Domain_Indices().Thickened(ghost_cells);dimension_start=domain_indices.Minimum_Corner();dimension_end=domain_indices.Maximum_Corner();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_UNIFORM<TV,T2>::
~EXTRAPOLATION_UNIFORM()
{} 
//#####################################################################
// Function Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_UNIFORM<TV,T2>::
Extrapolate(const T time,const bool fill_ghost_cells)
{
    ARRAY<T,TV_INT> phi_ghost(node_grid.Domain_Indices(ghost_cells),no_init);
    ARRAY<T2,TV_INT> u_ghost(node_grid.Domain_Indices(ghost_cells),no_init); 
    if(fill_ghost_cells){
        BOUNDARY<TV,T> phi_boundary;phi_boundary.Fill_Ghost_Cells(node_grid,phi,phi_ghost,0,time,ghost_cells);
        boundary->Fill_Ghost_Cells(node_grid,u,u_ghost,0,time,ghost_cells);}
    else{
        phi_ghost=phi;
        u_ghost=u;}
    if(isobaric_fix_width) phi_ghost+=isobaric_fix_width;   

    int heap_length=0;
    ARRAY<bool,TV_INT> done(node_grid.Domain_Indices(ghost_cells+1)); // extra cell for Update_Close_Point
    ARRAY<bool,TV_INT> close(node_grid.Domain_Indices(ghost_cells));
    ARRAY<TV_INT> heap(close.array.Size());
    heap.Fill(TV_INT()-1);
    Initialize(phi_ghost,done,close,heap,heap_length);

    while(heap_length && phi_ghost(heap(0)) <= band_width+isobaric_fix_width){
        TV_INT index=Remove_Root_From_Heap(phi_ghost,heap,heap_length,close);
        done(index)=true;
        Update_Close_Point(u_ghost,phi_ghost,done,index);
    
        if(collision_aware_extrapolation){
            for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && !close(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector))
                    Add_To_Heap(phi_ghost,heap,heap_length,close,index-axis_vector);
                if(index[axis] != dimension_end[axis]-1 && !done(index+axis_vector) && !close(index+axis_vector) && Neighbor_Visible(axis,index))
                    Add_To_Heap(phi_ghost,heap,heap_length,close,index+axis_vector);}}
        else{
            for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && !close(index-axis_vector))
                    Add_To_Heap(phi_ghost,heap,heap_length,close,index-axis_vector);
                if(index[axis] != dimension_end[axis]-1 && !done(index+axis_vector) && !close(index+axis_vector))
                    Add_To_Heap(phi_ghost,heap,heap_length,close,index+axis_vector);}}}

    ARRAYS_ND_BASE<T2,TV_INT>::Get(u,u_ghost);
    boundary->Apply_Boundary_Condition(node_grid,u,time);
}
//#####################################################################
// Function Initialize
//#####################################################################
// pass heap_length by reference
template<class TV,class T2> void EXTRAPOLATION_UNIFORM<TV,T2>::
Initialize(const T_ARRAYS_BASE& phi,ARRAYS_ND_BASE<bool,TV_INT>& done,ARRAYS_ND_BASE<bool,TV_INT>& close,ARRAY<TV_INT>& heap,int& heap_length)
{  
    assert(!seed_indices||!seed_done);
    if(seed_indices) for(int i=0;i<seed_indices->m;i++) done((*seed_indices)(i))=true;
    else if(seed_done) ARRAYS_ND_BASE<bool,TV_INT>::Put(*seed_done,done);
    else for(NODE_ITERATOR<TV> iterator(node_grid,ghost_cells);iterator.Valid();iterator.Next()) if(phi(iterator.Node_Index()) <= 0) done(iterator.Node_Index())=true;

    // find neighbors of done nodes which have positive phi
    if(collision_aware_extrapolation){
        for(NODE_ITERATOR<TV> iterator(node_grid,ghost_cells);iterator.Valid();iterator.Next()) if(done(iterator.Node_Index())){TV_INT index=iterator.Node_Index();
            for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && !close(index-axis_vector) && phi(index-axis_vector) > 0 && Neighbor_Visible(axis,index-axis_vector))
                    Add_To_Heap(phi,heap,heap_length,close,index-axis_vector);
                if(index[axis] != dimension_end[axis]-1 && !done(index+axis_vector) && !close(index+axis_vector) && phi(index+axis_vector) > 0 && Neighbor_Visible(axis,index))
                    Add_To_Heap(phi,heap,heap_length,close,index+axis_vector);}}}
    else{
        for(NODE_ITERATOR<TV> iterator(node_grid,ghost_cells);iterator.Valid();iterator.Next()) if(done(iterator.Node_Index())){TV_INT index=iterator.Node_Index();
            for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && !close(index-axis_vector) && phi(index-axis_vector) > 0)
                    Add_To_Heap(phi,heap,heap_length,close,index-axis_vector);
                if(index[axis] != dimension_end[axis]-1 && !done(index+axis_vector) && !close(index+axis_vector) && phi(index+axis_vector) > 0)
                    Add_To_Heap(phi,heap,heap_length,close,index+axis_vector);}}}
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
// needs done=0 around the outside of the domain
// note that sqrt(phix^2+phiy^2+phiz^2)=1 if it's a distance function
template<class TV,class T2> void EXTRAPOLATION_UNIFORM<TV,T2>::
Update_Close_Point(ARRAYS_ND_BASE<T2,TV_INT>& u,const T_ARRAYS_BASE& phi,const ARRAYS_ND_BASE<bool,TV_INT>& done,const TV_INT& index)
{
    T2 value[3]={}; // the value to use in the given direction
    T phix_dx[3]={0}; // the difference in phi value for the direction
    int number_of_axis=0; // the number of axis that we want to use later
    int missing_axis=2; // used in number_of_axis==2 case only, so it gives you which axis is missing (==2 for 2d)

    // check each principal axis
    for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis),low=index-axis_vector,high=index+axis_vector;
        bool check_low=done(low),check_high=done(high);
        if(collision_aware_extrapolation){
            if(check_low && !Neighbor_Visible(axis,low)) check_low=false;
            if(check_high && !Neighbor_Visible(axis,index)) check_high=false;}
        if(!check_low){
            if(check_high){value[number_of_axis]=u(high);phix_dx[number_of_axis]=phi(index)-phi(high);}
            else{missing_axis=axis;number_of_axis--;}}
        else if(!check_high){value[number_of_axis]=u(low);phix_dx[number_of_axis]=phi(index)-phi(low);}
        else{
            if(phi(low)<=phi(high)){value[number_of_axis]=u(low);phix_dx[number_of_axis]=phi(index)-phi(low);}
            else{value[number_of_axis]=u(high);phix_dx[number_of_axis]=phi(index)-phi(high);}}
        number_of_axis++;}

    assert(number_of_axis);
    if(number_of_axis==1) u(index)=value[0];
    else if(TV::m==2 || number_of_axis==2){
        T a=phix_dx[0]*optimization_scale[missing_axis],b=phix_dx[1],denominator=a+b,fraction=(T).5;
        if(denominator > small_number) fraction=clamp(a/denominator,(T)0,(T)1);
        u(index)=fraction*value[0]+(1-fraction)*value[1];}
    else{PHYSBAM_ASSERT(TV::m==3); // should only get here in 3D
        T a=phix_dx[0]*optimization_scale[1],b=phix_dx[1]*optimization_scale[0],c=phix_dx[2],denominator=a+b+c;
        T fraction_1=((T)1/3),fraction_2=((T)1/3);
        if(denominator > small_number){fraction_1=clamp(a/denominator,(T)0,(T)1);fraction_2=clamp(b/denominator,(T)0,1-fraction_1);}
        u(index)=fraction_1*value[0]+fraction_2*value[1]+((T)1-fraction_1-fraction_2)*value[2];}
}
//#####################################################################
namespace PhysBAM{
template class EXTRAPOLATION_UNIFORM<VECTOR<float,1>,float>;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,2>,float>;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,3>,float>;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,1>,VECTOR<float,1> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,1>,VECTOR<float,3> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,2>,VECTOR<float,2> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,2>,VECTOR<float,4> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,3>,VECTOR<float,3> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,3>,VECTOR<float,5> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,1>,double>;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,2>,double>;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,3>,double>;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,1>,VECTOR<double,1> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,1>,VECTOR<double,3> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,2>,VECTOR<double,2> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,2>,VECTOR<double,4> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,3>,VECTOR<double,3> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,3>,VECTOR<double,5> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2> >;
template class EXTRAPOLATION_UNIFORM<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3> >;
}
