//#####################################################################
// Copyright 2005-2008, Geoffrey Irving, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#ifdef USE_MPI
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_Tools/Parallel_Computation/LOCAL_GRID.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPI_UNIFORM_GRID<TV>::
MPI_UNIFORM_GRID(GRID<TV>& local_grid_input,const int number_of_ghost_cells_input,const bool skip_initialization,const TV_INT& processes_per_dimension,
    const VECTOR<bool,TV::m>& periodic_input,MPI::Group* group_input)
    :MPI_GRID<TV>(local_grid_input,number_of_ghost_cells_input,skip_initialization,processes_per_dimension,periodic_input,group_input)
{}
#ifdef USE_MPI
//#####################################################################
// Function Package_Cell_Data
//#####################################################################
template<class TV> template<class T2> MPI_PACKAGE MPI_UNIFORM_GRID<TV>::
Package_Cell_Data(ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& data,const RANGE<TV_INT>& region) const
{
    return MPI_PACKAGE(data,region);
}
//#####################################################################
// Function Package_Face_Data
//#####################################################################
template<class TV> template<class T_FACE_ARRAYS1> MPI_PACKAGE MPI_UNIFORM_GRID<TV>::
Package_Face_Data(T_FACE_ARRAYS1& data,const ARRAY<RANGE<TV_INT> >& regions) const
{
    MPI::Aint displacements[TV::m];MPI::Datatype old_types[TV::m];int lengths[TV::m];
    for(int axis=0;axis<TV::m;axis++){
        lengths[axis]=1;
        displacements[axis]=(MPI::Aint)&data.Component(axis)(regions(axis).Minimum_Corner());
        old_types[axis]=MPI_PACKAGE::Make_Arrays_Type(data.Component(axis),regions(axis));}
    MPI::Datatype datatype=MPI::Datatype::Create_struct(TV::m,lengths,displacements,old_types);
    for(int axis=0;axis<TV::m;axis++) old_types[axis].Free();
    return MPI_PACKAGE(datatype);
}
//#####################################################################
// Function Package_Common_Face_Data
//#####################################################################
template<class TV> template<class T_FACE_ARRAYS1> MPI_PACKAGE MPI_UNIFORM_GRID<TV>::
Package_Common_Face_Data(T_FACE_ARRAYS1& data,const int axis,const RANGE<TV_INT>& region) const
{
    return MPI_PACKAGE(data.Component(axis),region);
}
//#####################################################################
// Function Gather_Cell_Data
//#####################################################################
template<class TV> template<class T_ARRAYS> bool MPI_UNIFORM_GRID<TV>::
Gather_Cell_Data(const T_ARRAYS& local_data,T_ARRAYS& global_data) const
{
#if 0
    int tag=Get_Unique_Tag();
    int processes=comm->Get_size(),rank=comm->Get_rank(),master=0;
    GRID<TV> mac_global_grid=global_grid.Get_MAC_Grid();
    int ghost_cells=(local_grid.Domain_Indices().Minimum_Corner()-local_data.Domain_Indices().Minimum_Corner()).Max();

    RANGE<TV_INT> my_region=Find_Region_Box(rank+1,RANGE<TV_INT>::Zero_Box(),ghost_cells);
    GRID<TV> mac_local_grid=local_grid.Get_MAC_Grid();
    if(rank != master){
        ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
        MPI_PACKAGE package=Package_Cell_Data(const_cast<T_ARRAYS&>(local_data),my_region); // TODO change
        packages.Append(package);requests.Append(package.Isend(*comm,master,tag));
        MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);}
    else{
        LOCAL_GRID<TV> my_local_grid(mac_global_grid,mac_local_grid);
        my_local_grid.Put(local_data,my_region,global_data);
        for(int p=1;p<processes;p++){
            MPI::Status status;
            comm->Probe(MPI::ANY_SOURCE,tag,status);
            int source=status.Get_source();
            ARRAY<char> buffer(status.Get_count(MPI::PACKED));//int position=0;
            comm->Recv(&buffer(0),buffer.m,MPI::PACKED,source,tag);
            GRID<TV> other_grid=Restrict_Grid(all_coordinates(p));
            GRID<TV> mac_other_grid=other_grid.Get_MAC_Grid();
            T_ARRAYS other_array(mac_other_grid.Domain_Indices(ghost_cells));
            RANGE<TV_INT> other_region=Find_Region_Box(p,RANGE<TV_INT>::Zero_Box(),ghost_cells);
            MPI_PACKAGE package=Package_Cell_Data(other_array,other_region); // TODO: other_array.Domain_Indices());
            package.Unpack(buffer,*comm);
            package.Free();
            LOCAL_GRID<TV> other_local_grid(mac_global_grid,mac_other_grid);
            other_local_grid.Put(other_array,other_region,global_data);}}
    return rank!=master;
#endif
    return false;
}
//#####################################################################
// Function Scatter_Cell_Data
//#####################################################################
template<class TV> template<class T_ARRAYS> void MPI_UNIFORM_GRID<TV>::
Scatter_Cell_Data(const T_ARRAYS& global_data,T_ARRAYS& local_data) const
{
#if 0
    int tag=Get_Unique_Tag();
    int processes=comm->Get_size(),rank=comm->Get_rank(),master=0;
    int ghost_cells=(local_grid.Domain_Indices().Minimum_Corner()-local_data.Domain_Indices().Minimum_Corner()).Max();
    GRID<TV> mac_global_grid=global_grid.Get_MAC_Grid();
    GRID<TV> mac_local_grid=local_grid.Get_MAC_Grid();
    if(rank != master){
        MPI::Status status;
        comm->Probe(0,tag,status);
        ARRAY<char> buffer(status.Get_count(MPI::PACKED));
        comm->Recv(&buffer(0),buffer.m,MPI::PACKED,0,tag);
        MPI_PACKAGE package=Package_Cell_Data(local_data,local_data.Domain_Indices());
        package.Unpack(buffer,*comm);
        package.Free();}
    else{
        LOCAL_GRID<TV> my_local_grid(mac_global_grid,mac_local_grid);
        my_local_grid.Get(global_data,local_data);
        ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;ARRAY<T_ARRAYS > other_arrays(processes);
        for(int p=1;p<processes;p++){
            GRID<TV> other_grid=Restrict_Grid(all_coordinates(p));GRID<TV> mac_other_grid=other_grid.Get_MAC_Grid();
            LOCAL_GRID<TV> other_local_grid(mac_global_grid,mac_other_grid);
            other_arrays(p).Resize(mac_other_grid.Domain_Indices(ghost_cells),false,false);
            other_local_grid.Get(global_data,other_arrays(p));
            MPI_PACKAGE package=Package_Cell_Data(other_arrays(p),other_arrays(p).Domain_Indices());
            packages.Append(package);requests.Append(package.Isend(*comm,p-1,tag));}
        MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);}
#endif
}
//#####################################################################
// Function Get_Non_Overlapping_Face_Grid
//#####################################################################
template<class TV> GRID<TV> MPI_UNIFORM_GRID<TV>::
Get_Non_Overlapping_Face_Grid(const int axis) const
{
    GRID<TV> face_grid=local_grid.Get_Face_Grid(axis);
    if(side_neighbor_ranks(2*axis)!=MPI::PROC_NULL){
        const TV_INT counts=face_grid.Numbers_Of_Cells()+TV_INT::All_Ones_Vector()-TV_INT::Axis_Vector(axis);
        const RANGE<TV> box=face_grid.Domain()+RANGE<TV>(TV(),-face_grid.dX*TV::Axis_Vector(axis));
        return GRID<TV>(counts,box).Get_MAC_Grid_At_Regular_Positions();}
    else return face_grid.Get_MAC_Grid_At_Regular_Positions();
}
//#####################################################################

#else

//#####################################################################
namespace MPI{class Request{};}
namespace PhysBAM{class MPI_PACKAGE{};}
template<class TV> template<class T2> MPI_PACKAGE MPI_UNIFORM_GRID<TV>::Package_Cell_Data(ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >&,const RANGE<TV_INT>&) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_FACE_ARRAYS1> MPI_PACKAGE MPI_UNIFORM_GRID<TV>::Package_Face_Data(T_FACE_ARRAYS1&,const ARRAY<RANGE<TV_INT> >&) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_FACE_ARRAYS1> MPI_PACKAGE MPI_UNIFORM_GRID<TV>::Package_Common_Face_Data(T_FACE_ARRAYS1&,const int,const RANGE<TV_INT>&) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> GRID<TV> MPI_UNIFORM_GRID<TV>::Get_Non_Overlapping_Face_Grid(const int) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_ARRAYS> bool MPI_UNIFORM_GRID<TV>::Gather_Cell_Data(const T_ARRAYS& local_data,T_ARRAYS& global_data) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_ARRAYS> void MPI_UNIFORM_GRID<TV>::Scatter_Cell_Data(const T_ARRAYS& global_data,T_ARRAYS& local_data) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

//#####################################################################

#endif

//#####################################################################
namespace PhysBAM{
template GRID<VECTOR<double,1> > MPI_UNIFORM_GRID<VECTOR<double,1> >::Get_Non_Overlapping_Face_Grid(int) const;
template MPI_UNIFORM_GRID<VECTOR<double,1> >::MPI_UNIFORM_GRID(GRID<VECTOR<double,1> >&,int,bool,VECTOR<int,1> const&,VECTOR<bool,1> const&,MPI::Group*);
template GRID<VECTOR<double,2> > MPI_UNIFORM_GRID<VECTOR<double,2> >::Get_Non_Overlapping_Face_Grid(int) const;
template MPI_UNIFORM_GRID<VECTOR<double,2> >::MPI_UNIFORM_GRID(GRID<VECTOR<double,2> >&,int,bool,VECTOR<int,2> const&,VECTOR<bool,2> const&,MPI::Group*);
template GRID<VECTOR<double,3> > MPI_UNIFORM_GRID<VECTOR<double,3> >::Get_Non_Overlapping_Face_Grid(int) const;
template MPI_UNIFORM_GRID<VECTOR<double,3> >::MPI_UNIFORM_GRID(GRID<VECTOR<double,3> >&,int,bool,VECTOR<int,3> const&,VECTOR<bool,3> const&,MPI::Group*);
template GRID<VECTOR<float,1> > MPI_UNIFORM_GRID<VECTOR<float,1> >::Get_Non_Overlapping_Face_Grid(int) const;
template MPI_UNIFORM_GRID<VECTOR<float,1> >::MPI_UNIFORM_GRID(GRID<VECTOR<float,1> >&,int,bool,VECTOR<int,1> const&,VECTOR<bool,1> const&,MPI::Group*);
template GRID<VECTOR<float,2> > MPI_UNIFORM_GRID<VECTOR<float,2> >::Get_Non_Overlapping_Face_Grid(int) const;
template MPI_UNIFORM_GRID<VECTOR<float,2> >::MPI_UNIFORM_GRID(GRID<VECTOR<float,2> >&,int,bool,VECTOR<int,2> const&,VECTOR<bool,2> const&,MPI::Group*);
template GRID<VECTOR<float,3> > MPI_UNIFORM_GRID<VECTOR<float,3> >::Get_Non_Overlapping_Face_Grid(int) const;
template MPI_UNIFORM_GRID<VECTOR<float,3> >::MPI_UNIFORM_GRID(GRID<VECTOR<float,3> >&,int,bool,VECTOR<int,3> const&,VECTOR<bool,3> const&,MPI::Group*);
}
