//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Level_Sets/FAST_MARCHING.h>
#include <Geometry/Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FAST_MARCHING_METHOD_UNIFORM<TV>::
FAST_MARCHING_METHOD_UNIFORM(const LEVELSET<TV>& levelset,const int ghost_cells_input)
    :levelset(levelset),ghost_cells(ghost_cells_input),Neighbor_Visible(0)
{
    cell_grid=levelset.grid.Is_MAC_Grid()?levelset.grid:levelset.grid.Get_MAC_Grid_At_Regular_Positions();
    RANGE<TV_INT> domain_indices=cell_grid.Domain_Indices().Thickened(ghost_cells);dimension_start=domain_indices.Minimum_Corner();dimension_end=domain_indices.Maximum_Corner();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FAST_MARCHING_METHOD_UNIFORM<TV>::
~FAST_MARCHING_METHOD_UNIFORM()
{}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class TV> void FAST_MARCHING_METHOD_UNIFORM<TV>::
Fast_Marching_Method(ARRAY<T,TV_INT>& phi_ghost,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells,int process_sign)
{
    const RANGE<TV_INT>& domain=cell_grid.Domain_Indices(ghost_cells);
    int heap_length=0;
    ARRAY<T,TV_INT> phi_new(domain);
    for(CELL_ITERATOR<TV> iterator(cell_grid,domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();phi_new(cell)=phi_ghost(cell);}
    ARRAY<bool,TV_INT> done(domain.Thickened(1));
    ARRAY<int,TV_INT> close_k(domain.Thickened(1),false); // extra cell so that it is the same size as the done array for optimizations
    close_k.Fill(-1);
    ARRAY<TV_INT> heap(domain.Size(),false); // a generic version of heap((m+6)*(n+6)*(mn+6),false);
    Initialize_Interface(phi_new,done,close_k,heap,heap_length,seed_indices,add_seed_indices_for_ghost_cells);

    while(heap_length != 0){
        TV_INT index=heap(0); // smallest point is on top of heap
        if(stopping_distance && abs(phi_new(index)) > stopping_distance){ // exit early
            for(CELL_ITERATOR<TV> iterator(cell_grid,domain);iterator.Valid();iterator.Next()) if(!done(iterator.Cell_Index())){
                phi_new(iterator.Cell_Index())=LEVELSET_UTILITIES<T>::Sign(phi_new(iterator.Cell_Index()))*stopping_distance;}
            break;}
        done(index)=true;close_k(index)=-1; // add to done, remove from close
        FAST_MARCHING<T>::Down_Heap(phi_new,close_k,heap,heap_length);heap_length--; // remove point from heap

        if((process_sign==1 && phi_new(index)<0) || (process_sign==-1 && phi_new(index)>0)) continue;

        if(Neighbor_Visible)
            for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != domain.min_corner[axis] && !done(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector))
                    Update_Or_Add_Neighbor(phi_new,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != domain.max_corner[axis]-1 && !done(index+axis_vector) && Neighbor_Visible(axis,index))
                    Update_Or_Add_Neighbor(phi_new,done,close_k,heap,heap_length,index+axis_vector);}
        else
            for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != domain.min_corner[axis] && !done(index-axis_vector))
                    Update_Or_Add_Neighbor(phi_new,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != domain.max_corner[axis]-1 && !done(index+axis_vector))
                    Update_Or_Add_Neighbor(phi_new,done,close_k,heap,heap_length,index+axis_vector);}}

    RANGE<TV_INT> interior_domain=domain.Thickened(-ghost_cells);
    for(CELL_ITERATOR<TV> iterator(cell_grid,interior_domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();phi_ghost(cell)=phi_new(cell);}
}
//#####################################################################
// Function Update_Or_Add_Neighbor
//#####################################################################
template<class TV> inline void FAST_MARCHING_METHOD_UNIFORM<TV>::
Update_Or_Add_Neighbor(ARRAY<T,TV_INT>& phi_ghost,ARRAY<bool,TV_INT>& done,ARRAY<int,TV_INT>& close_k,ARRAY<TV_INT>& heap,int& heap_length,const TV_INT& neighbor)
{
    if(close_k(neighbor)>=0){
        Update_Close_Point(phi_ghost,done,neighbor);
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,close_k(neighbor));}
    else{
        close_k(neighbor)=0; // add to close 
        Update_Close_Point(phi_ghost,done,neighbor);
        heap(heap_length++)=neighbor;
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,heap_length-1);}
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
// pass heap_length by reference
template<class TV> void FAST_MARCHING_METHOD_UNIFORM<TV>::
Initialize_Interface(ARRAY<T,TV_INT>& phi_ghost,ARRAY<bool,TV_INT>& done,ARRAY<int,TV_INT>& close_k,ARRAY<TV_INT>& heap,int& heap_length,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells)
{ 
    if(seed_indices){
        for(int i=0;i<seed_indices->m;i++)Add_To_Initial(done,close_k,(*seed_indices)(i));
        if(add_seed_indices_for_ghost_cells){RANGE<TV_INT> ghost_domain=cell_grid.Domain_Indices().Thickened(ghost_cells);
            for(CELL_ITERATOR<TV> iterator(cell_grid,ghost_cells,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
                for(int i=0;i<GRID<TV>::number_of_neighbors_per_cell;i++){TV_INT neighbor_index(iterator.Cell_Neighbor(i));
                    if(ghost_domain.Lazy_Inside_Half_Open(neighbor_index) && LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(neighbor_index))){
                        if(!done(index))Add_To_Initial(done,close_k,index);
                        if(!done(neighbor_index))Add_To_Initial(done,close_k,neighbor_index);}}}}}
    else{
        ARRAY<T,TV_INT> phi_new(cell_grid.Domain_Indices(ghost_cells+1),false); // same size as done and close_k for array accelerations
        phi_new.Fill(2*cell_grid.dX.Max()); // ok positive, since minmag is used below
        if(Neighbor_Visible){
            for(FACE_ITERATOR<TV> iterator(cell_grid,ghost_cells,GRID<TV>::INTERIOR_REGION);iterator.Valid();iterator.Next()){TV_INT index1=iterator.First_Cell_Index(),index2=iterator.Second_Cell_Index();
                if(!Neighbor_Visible(iterator.Axis(),index1)){
                    if(phi_ghost(index1)<=0) Add_To_Initial(done,close_k,index1);
                    if(phi_ghost(index2)<=0) Add_To_Initial(done,close_k,index2);}
                else if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(index1),phi_ghost(index2))){
                    Add_To_Initial(done,close_k,index1);Add_To_Initial(done,close_k,index2);}}}
        else{
            for(FACE_ITERATOR<TV> iterator(cell_grid,ghost_cells,GRID<TV>::INTERIOR_REGION);iterator.Valid();iterator.Next()){TV_INT index1=iterator.First_Cell_Index(),index2=iterator.Second_Cell_Index();
                if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(index1),phi_ghost(index2))){
                    Add_To_Initial(done,close_k,index1);Add_To_Initial(done,close_k,index2);}}}

        LEVELSET<TV> levelset_ghost(cell_grid,phi_ghost);

        T fmm_initialization_iterative_tolerance_absolute=levelset.fmm_initialization_iterative_tolerance*cell_grid.dX.Min();    
 
        for(CELL_ITERATOR<TV> iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Cell_Index();
            if(!done(index)) continue;
            T value[3]={0}; // the phi value to use in the given direction
            int number_of_axis=0; // the number of axis that we want to use later
            TV location=iterator.Location();
            for(int axis=0;axis<TV::m;axis++){
                TV_INT axis_vector=TV_INT::Axis_Vector(axis),low=index-axis_vector,high=index+axis_vector;
                T dx=cell_grid.dX[axis];
                bool use_low=(done(low) && LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(low))),use_high=(done(high) && LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(high)));
                if(!use_low){
                    if(use_high) value[number_of_axis]=LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(high))*dx;
                    else number_of_axis--;}
                else if(!use_high) value[number_of_axis]=LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(low))*dx;
                else value[number_of_axis]=min(LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(high)),LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(low)))*dx;
                number_of_axis++;}
            if(number_of_axis==1) phi_new(index)=value[0];
            else if(number_of_axis==2){
                if(T d2=sqr(value[0])+sqr(value[1])) phi_new(index)=value[0]*value[1]/sqrt(d2);
                else phi_new(index)=0;}
            else{PHYSBAM_ASSERT(TV::m==3); // 2d should never get to this point
                T value_xy=value[0]*value[1],value_xz=value[0]*value[2],value_yz=value[1]*value[2],d2=sqr(value_xy)+sqr(value_xz)+sqr(value_yz);
                if(d2) phi_new(index)=value_xy*value[2]/sqrt(d2);
                else phi_new(index)=min(value[0],value[1],value[2]);}
            if(levelset.refine_fmm_initialization_with_iterative_solver){TV location=iterator.Location();
                TV vertex=location;T phi_value=phi_ghost(index);
                for(int iterations=0;iterations<levelset.fmm_initialization_iterations;iterations++){
                    if(abs(phi_value)>10*cell_grid.dX.Min())break; // stop if it looks like it's blowing up
                    vertex-=phi_value*levelset_ghost.Normal(vertex);phi_value=levelset_ghost.Phi(vertex);
                    if(abs(phi_value)<=fmm_initialization_iterative_tolerance_absolute){
                        T new_phi_value=(vertex-location).Magnitude();
                        if((new_phi_value-phi_new(index))/max(fmm_initialization_iterative_tolerance_absolute,phi_new(index))<levelset.fmm_initialization_iterative_drift_fraction && new_phi_value>0)
                            phi_new(index)=new_phi_value;
                        break;}}}
            phi_new(index)*=LEVELSET_UTILITIES<T>::Sign(phi_ghost(index));}

        for(CELL_ITERATOR<TV> iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(done(iterator.Cell_Index())) phi_ghost(iterator.Cell_Index())=phi_new(iterator.Cell_Index());} // initialize done points

    // initialize close points
    for(CELL_ITERATOR<TV> iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(close_k(iterator.Cell_Index())>=0){
        Update_Close_Point(phi_ghost,done,iterator.Cell_Index());
        heap(heap_length++)=iterator.Cell_Index();
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,heap_length-1);}
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
// pass heap_length by reference
template<class TV> void FAST_MARCHING_METHOD_UNIFORM<TV>::
Initialize_Interface(ARRAY<T,TV_INT>& phi_ghost,ARRAY<bool,TV_INT>& done,ARRAY<int,TV_INT>& close_k,ARRAY<TV_INT>& heap,int& heap_length,const bool add_seed_indices_for_ghost_cells)
{
    LEVELSET<TV> levelset_ghost(cell_grid,phi_ghost);

    for(CELL_ITERATOR<TV> iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(done(iterator.Cell_Index())) Add_To_Initial(done,close_k,iterator.Cell_Index());
    if(add_seed_indices_for_ghost_cells){RANGE<TV_INT> ghost_domain=cell_grid.Domain_Indices().Thickened(ghost_cells);
        for(CELL_ITERATOR<TV> iterator(cell_grid,ghost_cells,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
            for(int i=0;i<GRID<TV>::number_of_neighbors_per_cell;i++){TV_INT neighbor_index(iterator.Cell_Neighbor(i));
                if(ghost_domain.Lazy_Inside_Half_Open(neighbor_index) && LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(neighbor_index))){
                    if(!done(index))Add_To_Initial(done,close_k,index);
                    if(!done(neighbor_index))Add_To_Initial(done,close_k,neighbor_index);}}}}

   // initialize close points
    for(CELL_ITERATOR<TV> iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(close_k(iterator.Cell_Index())>=0){
        Update_Close_Point(phi_ghost,done,iterator.Cell_Index());
        heap(heap_length++)=iterator.Cell_Index();
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,heap_length-1);}
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
// needs done=0 around the outside of the domain
template<class TV> void FAST_MARCHING_METHOD_UNIFORM<TV>::
Update_Close_Point(ARRAY<T,TV_INT>& phi_ghost,const ARRAY<bool,TV_INT>& done,const TV_INT& index)
{
    T value[3]={}; // the phi value to use in the given direction
    T dx[3]={0}; // the edge length in the given direction
    int number_of_axis=0; // the number of axis that we want to use later

    // check each principal axis
    for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis),low=index-axis_vector,high=index+axis_vector;
        bool check_low=done(low),check_high=done(high);
        if(Neighbor_Visible){
            if(check_low && !Neighbor_Visible(axis,low)) check_low=false;
            if(check_high && !Neighbor_Visible(axis,index)) check_high=false;}
        dx[number_of_axis]=cell_grid.dX[axis];
        if(!check_low){
            if(check_high)value[number_of_axis]=phi_ghost(high);
            else number_of_axis--;}
        else if(!check_high) value[number_of_axis]=phi_ghost(low);
        else value[number_of_axis]=minmag(phi_ghost(low),phi_ghost(high));
        number_of_axis++;}

    phi_ghost(index)=FAST_MARCHING<T>::template Solve_Close_Point<TV::m>(phi_ghost(index),number_of_axis,value,dx);
}
//#####################################################################
// Function Add_To_Initial
//##################################################################### 
template<class TV> void FAST_MARCHING_METHOD_UNIFORM<TV>::
Add_To_Initial(ARRAY<bool,TV_INT>& done,ARRAY<int,TV_INT>& close_k,const TV_INT& index)
{
    done(index)=true;close_k(index)=-1; // add to done, remove from close 
    // add neighbors to close if not done 
    if(Neighbor_Visible){
        for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            if(!done(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector)) close_k(index-axis_vector)=0;
            if(!done(index+axis_vector) && Neighbor_Visible(axis,index)) close_k(index+axis_vector)=0;}}
    else{
        for(int axis=0;axis<TV::m;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            if(!done(index-axis_vector)) close_k(index-axis_vector)=0;
            if(!done(index+axis_vector)) close_k(index+axis_vector)=0;}}
}
//#####################################################################
namespace PhysBAM{
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<float,1> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<float,2> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<float,3> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<double,1> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<double,2> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<double,3> >;
}
