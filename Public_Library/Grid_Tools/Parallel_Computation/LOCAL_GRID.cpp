//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Parallel_Computation/LOCAL_GRID.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LOCAL_GRID<TV>::
LOCAL_GRID(const GRID<TV>& global_grid_input)
    :mpi_grid(grid,3,true),global_grid(global_grid_input)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LOCAL_GRID<TV>::
LOCAL_GRID(const GRID<TV>& global_grid_input,const GRID<TV>& local_grid_input)
    :grid(local_grid_input),mpi_grid(grid,3,true),global_grid(global_grid_input)
{
    Initialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LOCAL_GRID<TV>::
~LOCAL_GRID()
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void LOCAL_GRID<TV>::
Initialize()
{// find offset
    TV X=grid.Domain().Minimum_Corner();
    TV_INT local_offset=grid.Closest_Node(X),global_offset=global_grid.Closest_Node(X);
    offset=global_offset-local_offset;
    PHYSBAM_ASSERT((grid.Node(local_offset)-global_grid.Node(global_offset)).Magnitude()<(T).01*grid.dX.Min(),"mismatch between global and local grids");
    // pretend there are neighbors on all sides for use in Find_Boundary_Regions
    mpi_grid.side_neighbor_ranks.Resize(GRID<TV>::number_of_neighbors_per_cell);
    mpi_grid.side_neighbor_ranks.Fill(1);
    mpi_grid.all_neighbor_ranks.Resize(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    mpi_grid.all_neighbor_ranks.Fill(1);
    ARRAY<RANGE<TV_INT> > regions;
    mpi_grid.Find_Boundary_Regions(regions,RANGE<TV_INT>::Zero_Box(),true,RANGE<VECTOR<int,1> >(VECTOR<int,1>(-5),VECTOR<int,1>(-5)),true);
    RANGE<TV_INT> global_region=global_grid.Domain_Indices();
    neighbor_overlaps.Resize(regions.m);
    for(int n=0;n<regions.m;n++) if(regions(n).Lazy_Intersection(global_region-offset)) neighbor_overlaps(n)=true;
}
//#####################################################################
// Function Put
//#####################################################################
template<class TV> template<class T_ARRAYS> void LOCAL_GRID<TV>::
Put(const T_ARRAYS& local_data,const RANGE<TV_INT>& region,T_ARRAYS& global_data) const
{
    CELL_ITERATOR<TV> local(grid,region),global(grid,region+offset);
    for(;local.Valid();local.Next(),global.Next()) if(global_data.Valid_Index(global.Cell_Index()) && local_data.Valid_Index(local.Cell_Index()))
        global_data(global.Cell_Index())=local_data(local.Cell_Index());
}
//#####################################################################
// Function Put
//#####################################################################
template<class TV> template<class T_ARRAYS> void LOCAL_GRID<TV>::
Put(const T_ARRAYS& local_data,T_ARRAYS& global_data,const RANGE<TV_INT>& sentinels) const
{
    Put(local_data,Interior_Region(sentinels),global_data);
    ARRAY<RANGE<TV_INT> > regions;
    mpi_grid.Find_Boundary_Regions(regions,sentinels,true,RANGE<VECTOR<int,1> >(VECTOR<int,1>(-3),VECTOR<int,1>(-1)),true);
    for(int n=0;n<regions.m;n++) if(!neighbor_overlaps(n)) Put(local_data,regions(n),global_data);
}
//#####################################################################
// Function Get
//#####################################################################
template<class TV> template<class T_ARRAYS> void LOCAL_GRID<TV>::
Get(const T_ARRAYS& global_data,T_ARRAYS& local_data) const
{
    RANGE<TV_INT> region=local_data.Domain_Indices();
    CELL_ITERATOR<TV> local(grid,region),global(grid,region+offset);
    for(;local.Valid();local.Next(),global.Next()) local_data(local.Cell_Index())=global_data(global.Cell_Index());
}
//#####################################################################
// Function Put_Faces
//#####################################################################
template<class TV> template<class T_FACE_ARRAYS> void LOCAL_GRID<TV>::
Put_Faces(const T_FACE_ARRAYS& local_data,T_FACE_ARRAYS& global_data) const
{
    for(int axis=0;axis<TV::m;axis++) Put(local_data.Component(axis),global_data.Component(axis),mpi_grid.Face_Sentinels(axis));
}
//#####################################################################
// Function Maximum_Error
//#####################################################################
template<class TV> template<class T_ARRAYS> typename TV::SCALAR LOCAL_GRID<TV>::
Maximum_Error(const T_ARRAYS& local_data,const T_ARRAYS& global_data,const int bandwidth,TV_INT& index,const RANGE<TV_INT>& sentinels) const
{
    RANGE<TV_INT> region=Interior_Region(sentinels).Thickened(bandwidth);
    T max_error=0;CELL_ITERATOR<TV> local(grid,region),global(grid,region+offset);
    for(;local.Valid();local.Next(),global.Next()) if(global_data.Valid_Index(global.Cell_Index()) && local_data.Valid_Index(local.Cell_Index())){
            T error=sqrt((T)Magnitude_Squared(global_data(global.Cell_Index())-local_data(local.Cell_Index())));
        if(max_error<error){max_error=error;index=local.Cell_Index();}}
    return max_error;
}
//#####################################################################
// Function Maximum_Error
//#####################################################################
template<class TV> template<class T_FACE_ARRAYS> typename TV::SCALAR LOCAL_GRID<TV>::
Maximum_Error(const std::string& prefix,const T_FACE_ARRAYS& local_data,const T_FACE_ARRAYS& global_data,const int bandwidth,const T threshold)
{
    T max_error=(T)0;
    const char *axis_names[3]={"x","y","z"};
    for(int axis=0;axis<TV::m;axis++){
        TV_INT index;
        max_error=Maximum_Error(local_data.Component(axis),global_data.Component(axis),bandwidth,index,mpi_grid.Face_Sentinels(axis));
        if(max_error>threshold){LOG::cout<<prefix<<", face "<<axis_names[axis]<<" = "<<max_error<<" ("<<index<<")"<<std::endl;}
        }
    return max_error;
}
//#####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(T,d) \
    template class LOCAL_GRID<VECTOR<T,d> >; \
    template T LOCAL_GRID<VECTOR<T,d> >::Maximum_Error<ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > >(ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > const&,ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<VECTOR<T,d> >::Maximum_Error<ARRAY<VECTOR<T,d>,VECTOR<int,d> > >(ARRAY<VECTOR<T,d>,VECTOR<int,d> > const&,ARRAY<VECTOR<T,d>,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<VECTOR<T,d> >::Maximum_Error<ARRAY<T,VECTOR<int,d> > >(ARRAY<T,VECTOR<int,d> > const&,ARRAY<T,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<VECTOR<T,d> >::Maximum_Error<ARRAY<bool,VECTOR<int,d> > >(ARRAY<bool,VECTOR<int,d> > const&,ARRAY<bool,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<VECTOR<T,d> >::Maximum_Error<ARRAY<int,VECTOR<int,d> > >(ARRAY<int,VECTOR<int,d> > const&,ARRAY<int,VECTOR<int,d> > const&,int,VECTOR<int,d>&,RANGE<VECTOR<int,d> > const&) const; \
    template T LOCAL_GRID<VECTOR<T,d> >::Maximum_Error<ARRAY<bool,FACE_INDEX<d> > >(std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,ARRAY<bool,FACE_INDEX<d> > const&,ARRAY<bool,FACE_INDEX<d> > const&,int,T); \
    template T LOCAL_GRID<VECTOR<T,d> >::Maximum_Error<ARRAY<T,FACE_INDEX<d> > >(std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,ARRAY<T,FACE_INDEX<d> > const&,ARRAY<T,FACE_INDEX<d> > const&,int,T); \
    template void LOCAL_GRID<VECTOR<T,d> >::Put<ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > >(ARRAY<VECTOR<T,d+2>,VECTOR<int,d> > const&,ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<VECTOR<T,d> >::Put<ARRAY<VECTOR<T,d>,VECTOR<int,d> > >(ARRAY<VECTOR<T,d>,VECTOR<int,d> > const&,ARRAY<VECTOR<T,d>,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<VECTOR<T,d> >::Put<ARRAY<T,VECTOR<int,d> > >(ARRAY<T,VECTOR<int,d> > const&,ARRAY<T,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<VECTOR<T,d> >::Put<ARRAY<bool,VECTOR<int,d> > >(ARRAY<bool,VECTOR<int,d> > const&,ARRAY<bool,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<VECTOR<T,d> >::Put<ARRAY<int,VECTOR<int,d> > >(ARRAY<int,VECTOR<int,d> > const&,ARRAY<int,VECTOR<int,d> >&,RANGE<VECTOR<int,d> > const&) const; \
    template void LOCAL_GRID<VECTOR<T,d> >::Put_Faces<ARRAY<bool,FACE_INDEX<d> > >(ARRAY<bool,FACE_INDEX<d> > const&,ARRAY<bool,FACE_INDEX<d> >&) const; \
    template void LOCAL_GRID<VECTOR<T,d> >::Put_Faces<ARRAY<T,FACE_INDEX<d> > >(ARRAY<T,FACE_INDEX<d> > const&,ARRAY<T,FACE_INDEX<d> >&) const;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
}
