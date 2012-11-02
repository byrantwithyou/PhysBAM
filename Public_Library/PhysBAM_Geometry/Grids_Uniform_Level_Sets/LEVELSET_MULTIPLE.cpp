//#####################################################################
// Copyright 2005, Ron Fedkiw, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LEVELSET_MULTIPLE<T_GRID>::
LEVELSET_MULTIPLE(T_GRID& grid_input,ARRAY<T_ARRAYS_SCALAR>& phis_input,const bool use_external_levelsets_input)
    :grid(grid_input),phis(phis_input),levelset_callbacks(0),use_external_levelsets(use_external_levelsets_input)
{
    if(!use_external_levelsets) Recreate_Levelsets();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LEVELSET_MULTIPLE<T_GRID>::
~LEVELSET_MULTIPLE()
{
    if(!use_external_levelsets) for(int i=0;i<levelsets.m;i++) delete levelsets(i);
}
//#####################################################################
// Function Recreate_Levelsets
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Recreate_Levelsets() 
{
    assert(!use_external_levelsets);
    for(int i=0;i<levelsets.m;i++) delete levelsets(i);
    levelsets.Resize(phis.m);
    for(int i=0;i<levelsets.m;i++) levelsets(i)=new LEVELSET<TV>(grid,phis(i));
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Fill_Ghost_Cells(ARRAY<T_ARRAYS_SCALAR >& phi_ghost,const T time,const int number_of_ghost_cells) 
{
    for(int i=0;i<levelsets.m;i++) levelsets(i)->boundary->Fill_Ghost_Cells(grid,levelsets(i)->phi,phi_ghost(i),0,time,number_of_ghost_cells);
}
//#####################################################################
// Function Set_Custom_Boundary
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Set_Custom_Boundary(BOUNDARY<TV,T>& boundary_input)
{
    PHYSBAM_FATAL_ERROR(); // TODO: The next line does not compile.
//    for(int i=0;i<levelsets.m;i++) levelsets(i)->Set_Custom_Boundary(boundary_input);
}
//#####################################################################
// Function Set_Custom_Interpolation
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Set_Custom_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
{
    for(int i=0;i<levelsets.m;i++) levelsets(i)->Set_Custom_Interpolation(interpolation_input);
}
//#####################################################################
// Function Set_Custom_Secondary_Interpolation
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Set_Custom_Secondary_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
{
    for(int i=0;i<levelsets.m;i++) levelsets(i)->Set_Custom_Secondary_Interpolation(interpolation_input);
}
//#####################################################################
// Function Set_Custom_Normal_Interpolation
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Set_Custom_Normal_Interpolation(T_INTERPOLATION_VECTOR& normal_interpolation_input)
{
    for(int i=0;i<levelsets.m;i++) levelsets(i)->Set_Custom_Normal_Interpolation(normal_interpolation_input);
}
//#####################################################################
// Function Set_Custom_Curvature_Interpolation
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Set_Custom_Curvature_Interpolation(T_INTERPOLATION_SCALAR& curvature_interpolation_input)
{
    for(int i=0;i<levelsets.m;i++) levelsets(i)->Set_Custom_Curvature_Interpolation(curvature_interpolation_input);
}
//#####################################################################
// Function Set_Levelset_Callbacks
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Set_Levelset_Callbacks(LEVELSET_CALLBACKS<T_GRID>& levelset_callbacks_input)
{
    levelset_callbacks=&levelset_callbacks_input;
    for(int i=0;i<levelsets.m;i++) levelsets(i)->Set_Levelset_Callbacks(levelset_callbacks_input);
}
//#####################################################################
// Function Inside_Region
//#####################################################################
template<class T_GRID> int LEVELSET_MULTIPLE<T_GRID>::
Inside_Region(const TV_INT& index) const // assumes exactly one Phi<0 on a node
{
    for(int k=0;k<phis.m-1;k++) if(Phi(k,index)<=0) return k;
    assert(Phi(phis.m-1,index)<=0);
    return phis.m-1;
}
//#####################################################################
// Function Inside_Region
//#####################################################################
template<class T_GRID> int LEVELSET_MULTIPLE<T_GRID>::
Inside_Region(const TV_INT& index,T& phi) const // assumes exactly one Phi<0 on a node
{
    for(int k=0;k<phis.m-1;k++){phi=Phi(k,index);if(phi<=0) return k;}
    phi=Phi(phis.m-1,index);
    assert(phi<=0);
    return phis.m-1;
}
//#####################################################################
// Function Inside_Region
//#####################################################################
template<class T_GRID> int LEVELSET_MULTIPLE<T_GRID>::
Inside_Region(const TV& location) const
{
    T minimum_phi=Phi(0,location);
    int minimum_region=0;
    for(int k=1;k<phis.m;k++){T candidate_phi=Phi(k,location);if(candidate_phi<minimum_phi){minimum_phi=candidate_phi;minimum_region=k;}}
    return minimum_region;
}
//#####################################################################
// Function Inside_Region
//#####################################################################
template<class T_GRID> int LEVELSET_MULTIPLE<T_GRID>::
Inside_Region(const TV& location,T& phi) const
{
    T minimum_phi=Phi(0,location);
    int minimum_region=0;
    for(int k=1;k<phis.m;k++){T candidate_phi=Phi(k,location);if(candidate_phi<minimum_phi){minimum_phi=candidate_phi;minimum_region=k;}}
    phi=minimum_phi;
    return minimum_region;
}
//#####################################################################
// Function Inside_Region_Face
//#####################################################################
template<class T_GRID> int LEVELSET_MULTIPLE<T_GRID>::
Inside_Region_Face(const int axis,const TV_INT& face_index) const // does not assume exactly one Phi<0
{
    TV_INT cell_1,cell_2;
    grid.Cells_Touching_Face(axis,face_index,cell_1,cell_2);
    int region_1,region_2;
    T phi_1,phi_2;
    Minimum_Regions(cell_1,cell_2,region_1,region_2,phi_1,phi_2);
    if(phi_1<=phi_2) return region_1;
    return region_2;
}
//#####################################################################
// Function Two_Minimum_Regions
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Two_Minimum_Regions(const TV_INT& index,int& minimum_region,int& second_minimum_region,T& minimum_phi,T& second_minimum_phi) const
{
    T phi1=phis(0)(index),phi2=phis(1)(index);
    if(phi1<phi2){
        minimum_phi=phi1;
        minimum_region=0;
        second_minimum_phi=phi2;
        second_minimum_region=1;}
    else{
        minimum_phi=phi2;
        minimum_region=1;
        second_minimum_phi=phi1;
        second_minimum_region=0;}
    for(int k=2;k<phis.m;k++){
        T candidate_phi=phis(k)(index);
        if(candidate_phi<minimum_phi){
            second_minimum_phi=minimum_phi;
            second_minimum_region=minimum_region;
            minimum_phi=candidate_phi;
            minimum_region=k;}
        else if(candidate_phi<second_minimum_phi){
            second_minimum_phi=candidate_phi;
            second_minimum_region=k;}}
}
//#####################################################################
// Function Two_Minimum_Regions
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Two_Minimum_Regions(const TV& location,int& minimum_region,int& second_minimum_region,T& minimum_phi,T& second_minimum_phi) const
{
    T phi1=Phi(0,location),phi2=Phi(1,location);
    if(phi1<phi2){
        minimum_phi=phi1;
        minimum_region=0;
        second_minimum_phi=phi2;
        second_minimum_region=1;}
    else{
        minimum_phi=phi2;
        minimum_region=1;
        second_minimum_phi=phi1;
        second_minimum_region=0;}
    for(int k=2;k<phis.m;k++){
        T candidate_phi=Phi(k,location);
        if(candidate_phi<minimum_phi){
            second_minimum_phi=minimum_phi;
            second_minimum_region=minimum_region;
            minimum_phi=candidate_phi;
            minimum_region=k;}
        else if(candidate_phi<second_minimum_phi){
            second_minimum_phi=candidate_phi;
            second_minimum_region=k;}}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T::SCALAR LEVELSET_MULTIPLE<T_GRID>::
CFL(const T_FACE_ARRAYS_SCALAR& face_velocities) const
{
    T minimum_cfl=FLT_MAX;
    for(int i=0;i<levelsets.m;i++) minimum_cfl=min(minimum_cfl,levelsets(i)->CFL(face_velocities));
    return minimum_cfl;
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T::SCALAR LEVELSET_MULTIPLE<T_GRID>::
CFL(const ARRAY<TV,TV_INT>& velocity) const
{
    T minimum_cfl=FLT_MAX;
    for(int i=0;i<levelsets.m;i++) minimum_cfl=min(minimum_cfl,levelsets(i)->CFL(velocity));
    return minimum_cfl;
}
//#####################################################################
// Function Set_Collision_Body_List
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Set_Collision_Body_List(T_GRID_BASED_COLLISION_GEOMETRY& collision_body_list_input)
{
    for(int i=0;i<levelsets.m;i++) levelsets(i)->Set_Collision_Body_List(collision_body_list_input);
}
//#####################################################################
// Function Is_Projected_At_Nodes
//#####################################################################
template<class T_GRID> bool LEVELSET_MULTIPLE<T_GRID>::
Is_Projected_At_Nodes()
{
    for(T_CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        int inside=0;
        for(int i=0;i<levelsets.m;i++) if(Phi(i,iterator.Cell_Index())<=0) inside++;
        if(inside!=1) return false;}
    return true;
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Compute_Normals(const T time)
{
    for(int i=0;i<levelsets.m;i++) levelsets(i)->Compute_Normals(time);
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Compute_Curvature(const T time)
{
    for(int i=0;i<levelsets.m;i++) levelsets(i)->Compute_Curvature(time);
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Fast_Marching_Method(const ARRAY<int>& local_advection_spatial_orders,int process_sign)
{
    for(int i=0;i<levelsets.m;i++) levelsets(i)->Fast_Marching_Method(local_advection_spatial_orders(i),process_sign);
}
//#####################################################################
// Function Project_Levelset
//##################################################################### 
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Project_Levelset(const int number_of_ghost_cells)
{
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){
        int minimum_region,second_minimum_region;T minimum_phi,second_minimum_phi;
        Two_Minimum_Regions(iterator.Cell_Index(),minimum_region,second_minimum_region,minimum_phi,second_minimum_phi);
        T correction=(T).5*(minimum_phi+second_minimum_phi);
        for(int k=0;k<phis.m;k++) phis(k)(iterator.Cell_Index())-=correction;}
}
//#####################################################################
// Function Get_Single_Levelset
//#####################################################################
template<class T_GRID> void LEVELSET_MULTIPLE<T_GRID>::
Get_Single_Levelset(const ARRAY<bool>& positive_regions,LEVELSET<TV>& levelset,const bool flood_fill_for_bubbles)
{
    ARRAY<T,TV_INT>& phi_ghost=levelset.phi;phi_ghost.Resize(grid.Domain_Indices(3));
    if(flood_fill_for_bubbles){
        ARRAY<int,TV_INT> colors(grid.Domain_Indices(3));colors.Fill(-1);
        ARRAY<bool,FACE_INDEX<TV::m> > edge_is_blocked(grid,3);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
            if(!positive_regions(Inside_Region(iterator.Cell_Index()))) colors(iterator.Cell_Index())=-2;}
        FLOOD_FILL<TV::m> flood_fill;
        int number_of_colors=flood_fill.Flood_Fill(colors,edge_is_blocked);
        ARRAY<bool> color_touches_top_of_domain(number_of_colors);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(grid,0,T_GRID::BOUNDARY_REGION,4);iterator.Valid();iterator.Next()){
            if(colors(iterator.First_Cell_Index())>0)color_touches_top_of_domain(colors(iterator.First_Cell_Index()))=true;}
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
            T minimum_phi;Inside_Region(iterator.Cell_Index(),minimum_phi);
            if(colors(iterator.Cell_Index())>0 && color_touches_top_of_domain(colors(iterator.Cell_Index())))
                phi_ghost(iterator.Cell_Index())=-minimum_phi; // make levelset positive in the dirichlet regions
            else phi_ghost(iterator.Cell_Index())=minimum_phi;}}
    else{
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
            T minimum_phi;int minimum_region=Inside_Region(iterator.Cell_Index(),minimum_phi);
            if(positive_regions(minimum_region)) phi_ghost(iterator.Cell_Index())=-minimum_phi; // make levelset positive in the dirichlet regions
            else phi_ghost(iterator.Cell_Index())=minimum_phi;}}
    int number_of_positive_regions=0;for(int i=0;i<positive_regions.m;i++)if(positive_regions(i))number_of_positive_regions++;
    if(number_of_positive_regions>1)levelset.Fast_Marching_Method();
}
template class LEVELSET_MULTIPLE<GRID<VECTOR<float,1> > >;
template class LEVELSET_MULTIPLE<GRID<VECTOR<float,2> > >;
template class LEVELSET_MULTIPLE<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_MULTIPLE<GRID<VECTOR<double,1> > >;
template class LEVELSET_MULTIPLE<GRID<VECTOR<double,2> > >;
template class LEVELSET_MULTIPLE<GRID<VECTOR<double,3> > >;
#endif
