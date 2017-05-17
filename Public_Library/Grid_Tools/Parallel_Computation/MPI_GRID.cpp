//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_GRID
//#####################################################################
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Grid_Tools/Parallel_Computation/MPI_GRID.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#ifdef USE_MPI
#include <Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#endif
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <climits>
using namespace PhysBAM;

#ifdef USE_MPI

//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPI_GRID<TV>::
MPI_GRID(GRID<TV>& local_grid_input,const int number_of_ghost_cells_input,const bool skip_initialization,const TV_INT& processes_per_dimension,const VECTOR<bool,TV::m>& periodic_input,
    MPI::Group* group_input)
    :local_grid(local_grid_input),number_of_ghost_cells(number_of_ghost_cells_input),comm(0),group(group_input),current_tag(0),periodic(periodic_input),ignore_boundary_faces(false)
{
    if(skip_initialization) return;

    LOG::SCOPE scope("MPI INITIALIZE","Initializing MPI");
    if(!group) number_of_processes=MPI::COMM_WORLD.Get_size();
    else number_of_processes=group->Get_size();
    LOG::cout<<"number of processes = "<<number_of_processes<<std::endl;

    // extract global grid and divide among processes
    global_grid=local_grid.Get_Regular_Grid(); // Split_Grid and Initialize_Communicator currently assume grid is non-MAC
    Split_Grid(processes_per_dimension);

    // setup communicator and topology
    Initialize_Communicator(true,group);
    side_neighbor_ranks.Resize(GRID<TV>::number_of_neighbors_per_node);
    side_neighbor_directions.Resize(GRID<TV>::number_of_neighbors_per_node);
    all_neighbor_ranks.Resize(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    all_neighbor_directions.Resize(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    for(int n=0;n<GRID<TV>::number_of_neighbors_per_node;n++){
        side_neighbor_ranks(n)=process_ranks(GRID<TV>::Node_Neighbor(coordinates,n));
        side_neighbor_directions(n)=GRID<TV>::Node_Neighbor(TV_INT(),n);}
    for(int n=0;n<GRID<TV>::number_of_one_ring_neighbors_per_cell;n++){
        all_neighbor_ranks(n)=process_ranks(GRID<TV>::One_Ring_Neighbor(coordinates,n));
        all_neighbor_directions(n)=GRID<TV>::One_Ring_Neighbor(TV_INT(),n);}
    if(!group_input) group=new MPI::Group(comm->Get_group());

    // restrict this process to correct piece
    local_grid=Restrict_Grid(coordinates);
    // initialize offset
    TV_INT start_index;TV_INT end_index;
    for(int axis=0;axis<TV::m;axis++)
        start_index[axis]=boundaries(axis)(coordinates[axis]);
    local_to_global_offset=start_index-TV_INT::All_Ones_Vector();
    // initialize global column index boundaries
    ARRAY<VECTOR<int,2> > global_column_index_boundaries(number_of_processes);
    int offset=0;
    for(int proc=0;proc<all_coordinates.m;proc++){
        TV_INT proc_coordinates=all_coordinates(proc);
        for(int axis=0;axis<TV::m;axis++){
            start_index[axis]=boundaries(axis)(proc_coordinates[axis]);
            end_index[axis]=boundaries(axis)(proc_coordinates[axis]+1);}
        int owned_pressure_values=(end_index-start_index).Product();
        global_column_index_boundaries(proc)=VECTOR<int,2>(offset+1,offset+owned_pressure_values);
        offset+=owned_pressure_values;}
    // initialize the column indices of the ring of values that will be in the domain as well as one cell around it
    local_cell_index_to_global_column_index_map.Resize(local_grid.Domain_Indices(1));
    // first go through the interior and set what those n values will be
    offset=global_column_index_boundaries(rank+1).x;
    for(CELL_ITERATOR<TV> iterator(local_grid);iterator.Valid();iterator.Next()){
        local_cell_index_to_global_column_index_map(iterator.Cell_Index())=offset;
        offset++;}
    // now go through each of the boundaries and do those
    for(int axis=0;axis<TV::m;axis++)
        for(int axis_side=0;axis_side<2;axis_side++){
            int side=2*axis+axis_side;
            int neighbor_rank=side_neighbor_ranks(side);
            if(neighbor_rank>=0){
                TV_INT axis_vector=axis_side==0?-TV_INT::Axis_Vector(axis):TV_INT::Axis_Vector(axis);
                int start_column_index=global_column_index_boundaries(neighbor_rank+1).x;
                // Make that neighbors local_grid
                TV_INT proc_coordinates=all_coordinates(neighbor_rank+1);
                for(int temp_axis=0;temp_axis<TV::m;temp_axis++){
                    start_index[temp_axis]=boundaries(temp_axis)(proc_coordinates[temp_axis]);
                    end_index[temp_axis]=boundaries(temp_axis)(proc_coordinates[temp_axis]+1);}
                CELL_ITERATOR<TV> my_iterator(local_grid,0,GRID<TV>::BOUNDARY_INTERIOR_REGION,side);
                int neighbor_side=axis_side==0?2*axis+1:2*axis;
                GRID<TV> neighbor_grid=GRID<TV>(end_index-start_index+TV_INT::All_Ones_Vector(),RANGE<TV>(global_grid.X(start_index),global_grid.X(end_index))).Get_MAC_Grid();
                for(CELL_ITERATOR<TV> neighbor_iterator(neighbor_grid,0,GRID<TV>::BOUNDARY_INTERIOR_REGION,neighbor_side);neighbor_iterator.Valid();neighbor_iterator.Next(),my_iterator.Next()){
                    TV_INT my_cell_index=my_iterator.Cell_Index();
                    local_cell_index_to_global_column_index_map(my_cell_index+axis_vector)=neighbor_iterator.Flat_Index()+start_column_index-1;}}}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPI_GRID<TV>::
~MPI_GRID()
{
    if(comm){comm->Free();delete comm;}
    if(group){group->Free();delete group;}
}
//#####################################################################
// Function Initialize_Communicator
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Initialize_Communicator(const bool manual,MPI::Group* group)
{
    process_ranks.Resize(process_grid.Domain_Indices(1));process_ranks.array.Fill(MPI::PROC_NULL);
    TV_INT extents=process_grid.Domain_Indices().Maximum_Corner();
    comm=new MPI::Intracomm;
    if(!manual){ // setup cartesian communicator with the standard mpi function
        MPI::Cartcomm cartcomm=group?MPI::COMM_WORLD.Create(*group).Create_cart(TV::m,&extents[1],&periodic[1],true):
            MPI::COMM_WORLD.Create_cart(TV::m,&extents[1],&periodic[1],true);
        for(NODE_ITERATOR<TV> iterator(process_grid);iterator.Valid();iterator.Next()){TV_INT process=iterator.Node_Index();
            process_ranks(process)=cartcomm.Get_cart_rank(&(process-TV_INT::All_Ones_Vector())[1]);}
        *comm=cartcomm;}
    else{ // setup communicator manually to guarantee reasonable topology for SMP clusters
        // sort axes in decreasing order of how much we have to communicate along them
        ARRAY<int> axes(TV::m);ARRAY<T> axis_lengths(TV::m);
        for(int axis=0;axis<TV::m;axis++){axes(axis)=axis;axis_lengths(axis)=(T)global_grid.Domain_Indices().Maximum_Corner()[axis]/extents[axis];}
        axes.Sort([](int a,int b){return axis_lengths(a)<axis_lengths(b);});
        // lay out process ranks on grid
        Fill_Process_Ranks(process_grid,process_ranks,axes);
        // fill in ghost process_ranks for periodic domains
        if(periodic!=VECTOR<bool,TV::m>()) for(NODE_ITERATOR<TV> iterator(process_grid,1,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()){
            TV_INT node=iterator.Node_Index(),wrapped_node=node;
            for(int axis=0;axis<TV::m;axis++) if(periodic[axis]) wrapped_node[axis]=(node[axis]+process_grid.Counts()[axis])%process_grid.Counts()[axis];
            process_ranks(node)=process_ranks(wrapped_node);}
        // allocate communicator
        *comm=group?MPI::COMM_WORLD.Create(*group):MPI::COMM_WORLD.Dup();
        LOG::cout << "After Create comm " << std::endl;
    }
    all_coordinates.Resize(number_of_processes);
    for(NODE_ITERATOR<TV> iterator(process_grid);iterator.Valid();iterator.Next())
        all_coordinates(process_ranks(iterator.Node_Index())+1)=iterator.Node_Index();
    rank=comm->Get_rank();
    coordinates=all_coordinates(rank+1);
    LOG::cout<<"process_ranks = \n"<<process_ranks<<std::endl;
    LOG::cout<<"coordinates = "<<coordinates<<std::endl;
}
template<class T> static void Fill_Process_Ranks(GRID<VECTOR<T,1> >& process_grid,ARRAY<int,VECTOR<int,1> >& process_ranks,ARRAY<int>& axes)
{
    VECTOR<int,1> extents=process_grid.Domain_Indices().Maximum_Corner();
    for(int i=0;i<extents.x;i++)process_ranks(i)=i-1;
}
template<class T> static void Fill_Process_Ranks(GRID<VECTOR<T,2> >& process_grid,ARRAY<int,VECTOR<int,2> >& process_ranks,ARRAY<int>& axes)
{
    int next_rank=0;
    VECTOR<int,2> extents=process_grid.Domain_Indices().Maximum_Corner(),half_extents=extents/2;
    for(int i=0;i<half_extents[axes(0)];i++)for(int j=0;j<half_extents[axes(1)];j++)
        for(int ii=0;ii<2;ii++)for(int jj=0;jj<2;jj++){
            VECTOR<int,2> permuted_index(2*i+ii-1,2*j+jj-1),index;
            for(int a=0;a<2;a++)index[axes(a)]=permuted_index[a];
            process_ranks(index)=next_rank++;}
    for(int i=0;i<extents.x;i++)for(int j=0;j<extents.y;j++)if(process_ranks(i,j)==MPI::PROC_NULL) process_ranks(i,j)=next_rank++;
}
template<class T> static void Fill_Process_Ranks(GRID<VECTOR<T,3> >& process_grid,ARRAY<int,VECTOR<int,3> >& process_ranks,ARRAY<int>& axes)
{
    int next_rank=0;
    VECTOR<int,3> extents=process_grid.Domain_Indices().Maximum_Corner(),half_extents=extents/2;
    for(int i=0;i<half_extents[axes(0)];i++)for(int j=0;j<half_extents[axes(1)];j++)for(int ij=0;ij<half_extents[axes(2)];ij++)
        for(int ii=0;ii<2;ii++)for(int jj=0;jj<2;jj++)for(int ijij=0;ijij<2;ijij++){
            VECTOR<int,3> permuted_index(2*i+ii-1,2*j+jj-1,2*ij+ijij-1),index;
            for(int a=0;a<3;a++)index[axes(a)]=permuted_index[a];
            process_ranks(index)=next_rank++;}
    for(int i=0;i<extents.x;i++)for(int j=0;j<extents.y;j++)for(int ij=0;ij<extents.z;ij++)if(process_ranks(i,j,ij)==MPI::PROC_NULL) process_ranks(i,j,ij)=next_rank++;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Initialize(VECTOR<VECTOR<bool,2>,TV::m>& domain_walls)
{   
    // fix walls
    for(int i=0;i<GRID<TV>::number_of_neighbors_per_node;i++)
        if(side_neighbor_ranks(i)!=MPI::PROC_NULL) domain_walls(i/2)(i&1?0:1)=false;

    LOG::cout<<"mpi world rank = "<<MPI::COMM_WORLD.Get_rank()<<std::endl;
    LOG::cout<<"mpi cartesian rank = "<<rank<<std::endl;
    for(int axis=0;axis<TV::m;axis++)
        LOG::cout<<"mpi boundaries "<<axis<<" = "<<boundaries(axis)(coordinates[axis])<<" to "<<boundaries(axis)(coordinates[axis]+1)<<std::endl;
    LOG::cout<<"mpi topology = "<<process_grid.Domain_Indices()<<std::endl;
    LOG::cout<<"mpi process ranks = \n"<<process_ranks;
}
//#####################################################################
// Function Neighbor
//#####################################################################
template<class TV> bool MPI_GRID<TV>::
Neighbor(const int axis,const int axis_side) const
{   
    int side=2*axis+axis_side;
    return side_neighbor_ranks(side)!=MPI::PROC_NULL;
}
//#####################################################################

#endif

//#####################################################################
// Function Split_Grid
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Split_Grid(const TV_INT& processes_per_dimension)
{
    process_grid=GRID<TV>(Split_Grid(global_grid,processes_per_dimension),RANGE<TV>::Centered_Box());
}
//#####################################################################
// Function Split_Grid
//#####################################################################
template<class TV> VECTOR<int,1> MPI_GRID<TV>::
Split_Grid(const GRID<VECTOR<T,1> >& global_grid,const VECTOR<int,1>& processes_per_dimension)
{
    PHYSBAM_ASSERT(!processes_per_dimension.x || processes_per_dimension.x==number_of_processes);
    boundaries.Resize(1);
    Split_Dimension(global_grid.counts.x,number_of_processes,boundaries(0));
    return VECTOR<int,1>(number_of_processes);
}
//#####################################################################
// Function Split_Grid
//#####################################################################
template<class T>
static bool Minimize_2D_Surface_Area(const int number_of_processes,const int m,const int n,int& count_x)
{
    count_x=max(1,(int)sqrt((T)number_of_processes*m/n));
    if(number_of_processes%count_x==0){
        if(number_of_processes%(count_x+1)==0 && n*count_x+m*number_of_processes/count_x >= n*(count_x+1)+m*number_of_processes/(count_x+1)) count_x++;}
    else if(number_of_processes%(count_x+1)==0) count_x++;
    else return false;
    return true;
}
template<class TV> VECTOR<int,2> MPI_GRID<TV>::
Split_Grid(const GRID<VECTOR<T,2> >& global_grid,const VECTOR<int,2>& processes_per_dimension)
{
    int m=global_grid.counts.x,n=global_grid.counts.y;VECTOR<int,2> count;
    if(processes_per_dimension!=VECTOR<int,2>()) count=processes_per_dimension;
    else{ // try to figure out counts by minimizing surface area between processes
        if(!Minimize_2D_Surface_Area<T>(number_of_processes,m,n,count.x)){
            LOG::cerr<<"Don't know how to divide domain in both directions."<<std::endl;PHYSBAM_NOT_IMPLEMENTED();}
        count.y=number_of_processes/count.x;}
    PHYSBAM_ASSERT(count.x*count.y==number_of_processes);
    LOG::cout<<"dividing domain into "<<count.x<<" by "<<count.y<<" processor grid"<<std::endl;
    boundaries.Resize(1);
    Split_Dimension(m,count.x,boundaries(0));
    Split_Dimension(n,count.y,boundaries(1));
    return count;
}
//#####################################################################
// Function Split_Grid
//#####################################################################
template<class TV> VECTOR<int,3> MPI_GRID<TV>::
Split_Grid(const GRID<VECTOR<T,3> >& global_grid,const VECTOR<int,3>& processes_per_dimension)
{
    int m=global_grid.counts.x,n=global_grid.counts.y,mn=global_grid.counts.z;VECTOR<int,3> count;
    if(processes_per_dimension!=VECTOR<int,3>()) count=processes_per_dimension;
    else{ // try to figure out counts by minimizing surface area between processes
        int minimum_surface_area=INT_MAX;VECTOR<int,3> test_count;
        for(test_count.z=1;test_count.z<=number_of_processes;test_count.z++) if(number_of_processes%test_count.z==0){
            if(Minimize_2D_Surface_Area<T>(number_of_processes/test_count.z,m,n,test_count.x)){
                test_count.y=number_of_processes/(test_count.x*test_count.z);
                int surface_area=test_count.x*(n*mn)+test_count.y*(m*mn)+test_count.z*(m*n);
                if(surface_area<minimum_surface_area){count=test_count;minimum_surface_area=surface_area;}}
            if(Minimize_2D_Surface_Area<T>(number_of_processes/test_count.z,n,m,test_count.y)){
                test_count.x=number_of_processes/(test_count.y*test_count.z);
                int surface_area=test_count.x*(n*mn)+test_count.y*(m*mn)+test_count.z*(m*n);
                if(surface_area<minimum_surface_area){count=test_count;minimum_surface_area=surface_area;}}}
        if(minimum_surface_area==INT_MAX){LOG::cerr<<"Don't know how to divide domain in all directions."<<std::endl;PHYSBAM_NOT_IMPLEMENTED();}}
    PHYSBAM_ASSERT(count.x*count.y*count.z==number_of_processes);
    LOG::cout<<"dividing domain into "<<count.x<<" by "<<count.y<<" by "<<count.z<<" processor grid"<<std::endl;
    boundaries.Resize(3);
    Split_Dimension(m,count.x,boundaries(0));
    Split_Dimension(n,count.y,boundaries(1));
    Split_Dimension(mn,count.z,boundaries(2));
    return count;
}
//#####################################################################
// Function Split_Dimension
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Split_Dimension(const int m,const int processes,ARRAY<int>& boundaries)
{
    int cells_over_processes=m/processes,cells_mod_processes=m%processes;
    boundaries.Resize(processes+1);boundaries(0)=1;
    for(int p=0;p<processes;p++)boundaries(p+1)=boundaries(p)+cells_over_processes+(p<=cells_mod_processes);
}
//#####################################################################
// Function Restrict_Grid
//#####################################################################
template<class TV> GRID<TV> MPI_GRID<TV>::
Restrict_Grid(const TV_INT& coordinates) const
{
    TV_INT start_index,end_index;
    for(int axis=0;axis<TV::m;axis++){
        start_index[axis]=boundaries(axis)(coordinates[axis]);
        end_index[axis]=boundaries(axis)(coordinates[axis]+1);}
    return GRID<TV>(end_index-start_index+TV_INT::All_Ones_Vector(),RANGE<TV>(global_grid.X(start_index),global_grid.X(end_index))).Get_MAC_Grid();
}
//#####################################################################
// Function Find_Boundary_Regions
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,1> > >& regions,const RANGE<VECTOR<int,1> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
    const bool include_ghost_regions,const GRID<VECTOR<T,1> >& local_grid) const
{
    RANGE<VECTOR<int,1> > box=RANGE<VECTOR<int,1> >(VECTOR<int,1>(1),local_grid.numbers_of_cells)+sentinels;
    RANGE<VECTOR<int,1> > band_x=band;
    if(skip_common_boundary && band.max_corner.x>0){
        if(band.min_corner.x<0) PHYSBAM_NOT_IMPLEMENTED();
        if(sentinels.max_corner.x){band_x.min_corner.x++;band_x.max_corner.x++;}}
    regions.Clean_Memory();
    regions.Append(RANGE<VECTOR<int,1> >(box.min_corner,box.min_corner)+band_x);
    regions.Append(RANGE<VECTOR<int,1> >(box.max_corner,box.max_corner)-band_x);
}
//#####################################################################
// Function Find_Boundary_Regions
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,2> > >& regions,const RANGE<VECTOR<int,2> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
    const bool include_ghost_regions,const GRID<VECTOR<T,2> >& local_grid) const
{
    RANGE<VECTOR<int,2> > box=RANGE<VECTOR<int,2> >(VECTOR<int,2>(1,1),local_grid.numbers_of_cells)+sentinels;
    if(include_corners && include_ghost_regions){
        int band_size=band.Size();
        if(side_neighbor_ranks(0)<0) box.min_corner.x-=band_size;
        if(side_neighbor_ranks(1)<0) box.max_corner.x+=band_size;
        if(side_neighbor_ranks(2)<0) box.min_corner.y-=band_size;
        if(side_neighbor_ranks(3)<0) box.max_corner.y+=band_size;}
    RANGE<VECTOR<int,2> > band_x(VECTOR<int,2>(band.min_corner.x,0),VECTOR<int,2>(band.max_corner.x,0)),band_y(VECTOR<int,2>(0,band.min_corner.x),VECTOR<int,2>(0,band.max_corner.x));
    if(skip_common_boundary && band.max_corner.x>0){
        if(band.min_corner.x<0) PHYSBAM_NOT_IMPLEMENTED();
        if(sentinels.max_corner.x){band_x.min_corner.x++;band_x.max_corner.x++;}
        if(sentinels.max_corner.y){band_y.min_corner.y++;band_y.max_corner.y++;}}
    regions.Clean_Memory();
    if(!include_corners){ // add in side order
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.min_corner.x,box.min_corner.y),VECTOR<int,2>(box.min_corner.x,box.max_corner.y))+band_x);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.max_corner.x,box.min_corner.y),VECTOR<int,2>(box.max_corner.x,box.max_corner.y))-band_x);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.min_corner.x,box.min_corner.y),VECTOR<int,2>(box.max_corner.x,box.min_corner.y))+band_y);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.min_corner.x,box.max_corner.y),VECTOR<int,2>(box.max_corner.x,box.max_corner.y))-band_y);}
    else{ // add in one ring neighbors order
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.min_corner.x,box.min_corner.y),VECTOR<int,2>(box.min_corner.x,box.min_corner.y))+band_x+band_y);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.min_corner.x,box.min_corner.y),VECTOR<int,2>(box.max_corner.x,box.min_corner.y))+band_y);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.max_corner.x,box.min_corner.y),VECTOR<int,2>(box.max_corner.x,box.min_corner.y))-band_x+band_y);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.min_corner.x,box.min_corner.y),VECTOR<int,2>(box.min_corner.x,box.max_corner.y))+band_x);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.max_corner.x,box.min_corner.y),VECTOR<int,2>(box.max_corner.x,box.max_corner.y))-band_x);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.min_corner.x,box.max_corner.y),VECTOR<int,2>(box.min_corner.x,box.max_corner.y))+band_x-band_y);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.min_corner.x,box.max_corner.y),VECTOR<int,2>(box.max_corner.x,box.max_corner.y))-band_y);
        regions.Append(RANGE<VECTOR<int,2> >(VECTOR<int,2>(box.max_corner.x,box.max_corner.y),VECTOR<int,2>(box.max_corner.x,box.max_corner.y))-band_x-band_y);}
}
//#####################################################################
// Function Find_Boundary_Regions
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,3> > >& regions,const RANGE<VECTOR<int,3> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
    const bool include_ghost_regions,const GRID<VECTOR<T,3> >& local_grid) const
{
    RANGE<VECTOR<int,3> > box=RANGE<VECTOR<int,3> >(VECTOR<int,3>(1,1,1),local_grid.numbers_of_cells)+sentinels;
    if(include_corners && include_ghost_regions){
        int band_size=band.Size();
        if(side_neighbor_ranks(0)<0) box.min_corner.x-=band_size;
        if(side_neighbor_ranks(1)<0) box.max_corner.x+=band_size;
        if(side_neighbor_ranks(2)<0) box.min_corner.y-=band_size;
        if(side_neighbor_ranks(3)<0) box.max_corner.y+=band_size;
        if(side_neighbor_ranks(4)<0) box.min_corner.z-=band_size;
        if(side_neighbor_ranks(5)<0) box.max_corner.z+=band_size;}
    RANGE<VECTOR<int,3> > band_x(VECTOR<int,3>(band.min_corner.x,0,0),VECTOR<int,3>(band.max_corner.x,0,0)),band_y(VECTOR<int,3>(0,band.min_corner.x,0),VECTOR<int,3>(0,band.max_corner.x,0)),band_z(VECTOR<int,3>(0,0,band.min_corner.x),VECTOR<int,3>(0,0,band.max_corner.x));
    if(skip_common_boundary && band.max_corner.x>0){
        if(band.min_corner.x<0) PHYSBAM_NOT_IMPLEMENTED();
        if(sentinels.max_corner.x){band_x.min_corner.x++;band_x.max_corner.x++;}
        if(sentinels.max_corner.y){band_y.min_corner.y++;band_y.max_corner.y++;}
        if(sentinels.max_corner.z){band_z.min_corner.z++;band_z.max_corner.z++;}}
    regions.Clean_Memory();
    if(!include_corners){ // add in side order
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.max_corner.z))+band_x);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_x);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.max_corner.z))+band_y);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_y);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.min_corner.z))+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.max_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_z);}
    else{ // add in one ring neighbors order
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z))+band_x+band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.min_corner.z))+band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.min_corner.z))-band_x+band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.min_corner.z))+band_x+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.min_corner.z))+band_z); 
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.min_corner.z))-band_x+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.min_corner.z),VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.min_corner.z))+band_x-band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.min_corner.z))-band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.min_corner.z))-band_x-band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.max_corner.z))+band_x+band_y);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.max_corner.z))+band_y);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.max_corner.z))-band_x+band_y);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.max_corner.z))+band_x);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_x);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.min_corner.z),VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.max_corner.z))+band_x-band_y);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_y);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.min_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_x-band_y);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.max_corner.z),VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.max_corner.z))+band_x+band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.max_corner.z),VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.max_corner.z))+band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.max_corner.z),VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.max_corner.z))-band_x+band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.max_corner.z),VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.max_corner.z))+band_x-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.min_corner.y,box.max_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_z); 
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.min_corner.y,box.max_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_x-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.max_corner.z),VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.max_corner.z))+band_x-band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.min_corner.x,box.max_corner.y,box.max_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z),VECTOR<int,3>(box.max_corner.x,box.max_corner.y,box.max_corner.z))-band_x-band_y-band_z);
    }
}
//#####################################################################
// Function Find_Boundary_Regions
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Find_Boundary_Regions(ARRAY<RANGE<TV_INT> >& regions,const RANGE<TV_INT>& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
    const bool include_ghost_regions) const
{
    Find_Boundary_Regions(regions,sentinels,skip_common_boundary,band,include_corners,include_ghost_regions,local_grid);
}
//#####################################################################
// Function Find_Region_Box
//#####################################################################
template<class T> RANGE<VECTOR<int,1> >
Find_Region_Box_Helper(const ARRAY<int,VECTOR<int,1> >& process_ranks,const GRID<VECTOR<T,1> >& other_grid,const VECTOR<int,1>& coordinates,const RANGE<VECTOR<int,1> >& sentinels,const int band)
{
    RANGE<VECTOR<int,1> > box=RANGE<VECTOR<int,1> >(VECTOR<int,1>(),other_grid.numbers_of_cells)+sentinels;
    if(process_ranks(GRID<VECTOR<T,1> >::Node_Neighbor(coordinates,0))<0) box.min_corner.x-=band;
    if(process_ranks(GRID<VECTOR<T,1> >::Node_Neighbor(coordinates,1))<0) box.max_corner.x+=band;
    return box;
}
//#####################################################################
// Function Find_Region_Box
//#####################################################################
template<class T> RANGE<VECTOR<int,2> >
Find_Region_Box_Helper(const ARRAY<int,VECTOR<int,2> >& process_ranks,const GRID<VECTOR<T,2> >& other_grid,const VECTOR<int,2>& coordinates,const RANGE<VECTOR<int,2> >& sentinels,const int band)
{
    RANGE<VECTOR<int,2> > box=RANGE<VECTOR<int,2> >(VECTOR<int,2>(),other_grid.numbers_of_cells)+sentinels;
    if(process_ranks(GRID<VECTOR<T,2> >::Node_Neighbor(coordinates,0))<0) box.min_corner.x-=band;
    if(process_ranks(GRID<VECTOR<T,2> >::Node_Neighbor(coordinates,1))<0) box.max_corner.x+=band;
    if(process_ranks(GRID<VECTOR<T,2> >::Node_Neighbor(coordinates,2))<0) box.min_corner.y-=band;
    if(process_ranks(GRID<VECTOR<T,2> >::Node_Neighbor(coordinates,3))<0) box.max_corner.y+=band;
    return box;
}
//#####################################################################
// Function Find_Region_Box
//#####################################################################
template<class T> RANGE<VECTOR<int,3> >
Find_Region_Box_Helper(const ARRAY<int,VECTOR<int,3> >& process_ranks,const GRID<VECTOR<T,3> >& other_grid,const VECTOR<int,3>& coordinates,const RANGE<VECTOR<int,3> >& sentinels,const int band)
{
    RANGE<VECTOR<int,3> > box=RANGE<VECTOR<int,3> >(VECTOR<int,3>(),other_grid.numbers_of_cells)+sentinels;
    if(process_ranks(GRID<VECTOR<T,3> >::Node_Neighbor(coordinates,0))<0) box.min_corner.x-=band;
    if(process_ranks(GRID<VECTOR<T,3> >::Node_Neighbor(coordinates,1))<0) box.max_corner.x+=band;
    if(process_ranks(GRID<VECTOR<T,3> >::Node_Neighbor(coordinates,2))<0) box.min_corner.y-=band;
    if(process_ranks(GRID<VECTOR<T,3> >::Node_Neighbor(coordinates,3))<0) box.max_corner.y+=band;
    if(process_ranks(GRID<VECTOR<T,3> >::Node_Neighbor(coordinates,4))<0) box.min_corner.z-=band;
    if(process_ranks(GRID<VECTOR<T,3> >::Node_Neighbor(coordinates,5))<0) box.max_corner.z+=band;
    return box;
}
//#####################################################################
// Function Find_Region_Box
//#####################################################################
template<class TV> RANGE<VECTOR<int,TV::m>> MPI_GRID<TV>::
Find_Region_Box(const int processor,const RANGE<TV_INT>& sentinels,const int band) const
{
    TV_INT other_coordinates=all_coordinates(processor);
    GRID<TV> other_grid=Restrict_Grid(other_coordinates);
    return Find_Region_Box_Helper<T>(process_ranks,other_grid,other_coordinates,sentinels,band);
}
//#####################################################################

#ifdef USE_MPI

//#####################################################################
// Function Synchronize_Dt
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Synchronize_Dt(T& dt) const
{
    T dt_local=dt;
    MPI_UTILITIES::Reduce(dt_local,dt,MPI::MIN,*comm);
}
//#####################################################################
// Function Synchronize_J_Bounds
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Synchronize_J_Bounds(int& jmin,int& jmax) const
{
    int input_bounds[2]={jmin,-jmax},output_bounds[2];
    comm->Allreduce(input_bounds,output_bounds,2,MPI::INT,MPI::MIN);
    jmin=output_bounds[0];jmax=-output_bounds[1];
}
//#####################################################################
// Function Exchange_Boundary_Cell_Data
//#####################################################################
template<class TV> template<class T_MPI_GRID,class T2> void MPI_GRID<TV>::
Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& data,const int bandwidth,const bool include_corners) const
{
    PHYSBAM_ASSERT(bandwidth>0,"0 bandwidth exchange");
    Exchange_Boundary_Cell_Data(mpi_grid,local_grid,data,bandwidth,include_corners);
}
template<class TV> template<class T_MPI_GRID,class T2,class T_ARRAYS,class INDEX> void MPI_GRID<TV>::
Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,ARRAY_BASE<T2,T_ARRAYS,INDEX>& data,const int bandwidth,const bool include_corners) const
{
    PHYSBAM_ASSERT(bandwidth>0,"0 bandwidth exchange");
    Exchange_Boundary_Cell_Data(mpi_grid,local_grid,data,bandwidth,include_corners);
}
//#####################################################################
// Function Exchange_Boundary_Cell_Data
//#####################################################################
template<class TV> template<class T_MPI_GRID,class T2> void MPI_GRID<TV>::
Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,const GRID<TV>& local_grid,ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& data,const int bandwidth,const bool include_corners) const
{
    RANGE<TV_INT> sentinels=RANGE<TV_INT>::Zero_Box();
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    const ARRAY<int>& neighbor_ranks=include_corners?all_neighbor_ranks:side_neighbor_ranks;
    const ARRAY<TV_INT>& neighbor_directions=include_corners?all_neighbor_directions:side_neighbor_directions;
    // send
    ARRAY<RANGE<TV_INT> > send_regions;Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),include_corners,true,local_grid);
    for(int n=0;n<send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL){
        MPI_PACKAGE package=mpi_grid.Package_Cell_Data(data,send_regions(n));
        packages.Append(package);requests.Append(package.Isend(*comm,neighbor_ranks(n),Get_Send_Tag(neighbor_directions(n))));}
    // receive
    ARRAY<RANGE<TV_INT> > recv_regions;Find_Boundary_Regions(recv_regions,sentinels,false,RANGE<VECTOR<int,1> >(-bandwidth,-1),include_corners,true,local_grid);
    for(int n=0;n<recv_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL){
        MPI_PACKAGE package=mpi_grid.Package_Cell_Data(data,recv_regions(n));
        packages.Append(package);requests.Append(package.Irecv(*comm,neighbor_ranks(n),Get_Recv_Tag(neighbor_directions(n))));}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
}
template<class TV> template<class T_MPI_GRID,class T2,class T_ARRAYS,class INDEX> void MPI_GRID<TV>::
Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,const GRID<TV>& local_grid,ARRAY_BASE<T2,T_ARRAYS,INDEX>& data,const int bandwidth,const bool include_corners) const
{
    RANGE<TV_INT> sentinels=RANGE<TV_INT>::Zero_Box();
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    const ARRAY<int>& neighbor_ranks=include_corners?all_neighbor_ranks:side_neighbor_ranks;
    const ARRAY<TV_INT>& neighbor_directions=include_corners?all_neighbor_directions:side_neighbor_directions;
    // send
    ARRAY<RANGE<TV_INT> > send_regions;Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),include_corners,true,local_grid);
    for(int n=0;n<send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL){
        MPI_PACKAGE package=mpi_grid.Package_Cell_Data(data,send_regions(n));
        packages.Append(package);requests.Append(package.Isend(*comm,neighbor_ranks(n),Get_Send_Tag(neighbor_directions(n))));}
    // receive
    ARRAY<RANGE<TV_INT> > recv_regions;Find_Boundary_Regions(recv_regions,sentinels,false,RANGE<VECTOR<int,1> >(-bandwidth,-1),include_corners,true,local_grid);
    for(int n=0;n<recv_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL){
        MPI_PACKAGE package=mpi_grid.Package_Cell_Data(data,recv_regions(n));
        packages.Append(package);requests.Append(package.Irecv(*comm,neighbor_ranks(n),Get_Recv_Tag(neighbor_directions(n))));}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Exchange_Boundary_Face_Data
//#####################################################################
template<class TV> template<class T_MPI_GRID,class ARRAY<T,FACE_INDEX<TV::m> >> void MPI_GRID<TV>::
Exchange_Boundary_Face_Data(const T_MPI_GRID& mpi_grid,ARRAY<T,FACE_INDEX<TV::m> >& data,const int bandwidth) const
{   
    RANGE<VECTOR<int,1> > boundary_band(0,bandwidth-1),ghost_band(-bandwidth,-1);
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    // send
    ARRAY<ARRAY<RANGE<TV_INT> > > send_regions(T_MPI_GRID::GRID_T::dimension);
    for(int axis=0;axis<T_MPI_GRID::GRID_T::dimension;axis++)Find_Boundary_Regions(send_regions(axis),mpi_grid.Face_Sentinels(axis),true,boundary_band,true);
    for(int n=0;n<send_regions(0).m;n++)if(all_neighbor_ranks(n)!=MPI::PROC_NULL){
        ARRAY<RANGE<TV_INT> > send_regions_n(T_MPI_GRID::GRID_T::dimension);for(int axis=0;axis<T_MPI_GRID::GRID_T::dimension;axis++)send_regions_n(axis)=send_regions(axis)(n);
        MPI_PACKAGE package=mpi_grid.Package_Face_Data(data,send_regions_n);
        packages.Append(package);requests.Append(package.Isend(*comm,all_neighbor_ranks(n),Get_Send_Tag(all_neighbor_directions(n))));}
    // receive
    ARRAY<ARRAY<RANGE<TV_INT> > > recv_regions(T_MPI_GRID::GRID_T::dimension);
    for(int axis=0;axis<T_MPI_GRID::GRID_T::dimension;axis++)Find_Boundary_Regions(recv_regions(axis),mpi_grid.Face_Sentinels(axis),true,ghost_band,true);
    for(int n=0;n<recv_regions(0).m;n++)if(all_neighbor_ranks(n)!=MPI::PROC_NULL){
        ARRAY<RANGE<TV_INT> > recv_regions_n(T_MPI_GRID::GRID_T::dimension);for(int axis=0;axis<T_MPI_GRID::GRID_T::dimension;axis++)recv_regions_n(axis)=recv_regions(axis)(n);
        MPI_PACKAGE package=mpi_grid.Package_Face_Data(data,recv_regions_n);
        packages.Append(package);requests.Append(package.Irecv(*comm,all_neighbor_ranks(n),Get_Recv_Tag(all_neighbor_directions(n))));}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Average_Common_Face_Data
//#####################################################################
template<class TV> template<class T_MPI_GRID,class ARRAY<T,FACE_INDEX<TV::m> >> void MPI_GRID<TV>::
Average_Common_Face_Data(const T_MPI_GRID& mpi_grid,ARRAY<T,FACE_INDEX<TV::m> >& data) const
{
    //PHYSBAM_NOT_IMPLEMENTED(); // TODO: we probably don't need this function, but I'll leave it here for now

    ARRAY<ARRAY<RANGE<TV_INT> > > regions(TV::m);
    for(int axis=0;axis<TV::m;axis++)Find_Boundary_Regions(regions(axis),mpi_grid.Parallel_Face_Sentinels(axis),false,RANGE<VECTOR<int,1> >(0,0),false);
    ARRAY<MPI_PACKAGE> packages(regions(0).m);ARRAY<ARRAY<T> > buffers(regions(0).m);
    MPI::Datatype T_type=MPI_UTILITIES::Datatype<T>();
    // send and receive into temporary buffers
    ARRAY<MPI::Request> requests;
    for(int n=0;n<regions(0).m;n++)if(side_neighbor_ranks(n)!=MPI::PROC_NULL){int axis=n/2;
        packages(n)=mpi_grid.Package_Common_Face_Data(data,axis,regions(axis)(n));
        requests.Append(packages(n).Isend(*comm,side_neighbor_ranks(n),Get_Send_Tag(side_neighbor_directions(n))));
        buffers(n).Resize(packages(n).Size()/sizeof(T));
        requests.Append(comm->Irecv(buffers(n).Get_Array_Pointer(),buffers(n).m,T_type,side_neighbor_ranks(n),Get_Recv_Tag(side_neighbor_directions(n))));}
    // wait
    MPI_UTILITIES::Wait_All(requests);
    // average received data with local data (TODO: find a cleaner general way to do this)
    for(int n=0;n<regions(0).m;n++)if(side_neighbor_ranks(n)!=MPI::PROC_NULL){
        ARRAY<char> pack_buffer(packages(n).Pack_Size(*comm));packages(n).Pack(pack_buffer,*comm);
        ARRAY<T> local_buffer(buffers(n).m);
        int position=0;T_type.Unpack(pack_buffer.Get_Array_Pointer(),pack_buffer.m,local_buffer.Get_Array_Pointer(),local_buffer.m,position,*comm);
        local_buffer.Copy((T).5,buffers(n),(T).5,local_buffer); // average
        position=0;T_type.Pack(local_buffer.Get_Array_Pointer(),local_buffer.m,pack_buffer.Get_Array_Pointer(),pack_buffer.m,position,*comm);
        packages(n).Unpack(pack_buffer,*comm);}
    MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Union_Common_Face_Data
//#####################################################################
template<class TV> template<class T_MPI_GRID,class ARRAY<bool,FACE_INDEX<TV::m> > > void MPI_GRID<TV>::
Union_Common_Face_Data(const T_MPI_GRID& mpi_grid,ARRAY<bool,FACE_INDEX<TV::m> >& data) const
{
    ARRAY<ARRAY<RANGE<TV_INT> > > regions(TV::m);
    for(int axis=0;axis<TV::m;axis++)Find_Boundary_Regions(regions(axis),mpi_grid.Parallel_Face_Sentinels(axis),false,RANGE<VECTOR<int,1> >(0,0),false);
    ARRAY<MPI_PACKAGE> packages(regions(0).m);ARRAY<ARRAY<bool> > buffers(regions(0).m);
    MPI::Datatype T_type=MPI_UTILITIES::Datatype<bool>();
    // send and receive into temporary buffers
    ARRAY<MPI::Request> requests;
    for(int n=0;n<regions(0).m;n++)if(side_neighbor_ranks(n)!=MPI::PROC_NULL){int axis=n/2;
        packages(n)=mpi_grid.Package_Common_Face_Data(data,axis,regions(axis)(n));
        requests.Append(packages(n).Isend(*comm,side_neighbor_ranks(n),Get_Send_Tag(side_neighbor_directions(n))));
        buffers(n).Resize(packages(n).Size()/sizeof(bool));
        requests.Append(comm->Irecv(buffers(n).Get_Array_Pointer(),buffers(n).m,T_type,side_neighbor_ranks(n),Get_Recv_Tag(side_neighbor_directions(n))));}
    // wait
    MPI_UTILITIES::Wait_All(requests);
    // average received data with local data (TODO: find a cleaner general way to do this)
    for(int n=0;n<regions(0).m;n++)if(side_neighbor_ranks(n)!=MPI::PROC_NULL){
        ARRAY<char> pack_buffer(packages(n).Pack_Size(*comm));packages(n).Pack(pack_buffer,*comm);
        ARRAY<bool> local_buffer(buffers(n).m);
        int position=0;T_type.Unpack(pack_buffer.Get_Array_Pointer(),pack_buffer.m,local_buffer.Get_Array_Pointer(),local_buffer.m,position,*comm);
        for(int i=0;i<buffers(n).m;i++) local_buffer(i)=local_buffer(i)||buffers(n)(i);
        position=0;T_type.Pack(local_buffer.Get_Array_Pointer(),local_buffer.m,pack_buffer.Get_Array_Pointer(),pack_buffer.m,position,*comm);
        packages(n).Unpack(pack_buffer,*comm);}
    MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Copy_Common_Face_Data
//#####################################################################
template<class TV> template<class T_MPI_GRID,class ARRAY<T,FACE_INDEX<TV::m> >> void MPI_GRID<TV>::
Copy_Common_Face_Data(const T_MPI_GRID& mpi_grid,ARRAY<T,FACE_INDEX<TV::m> >& data) const
{
    ARRAY<ARRAY<RANGE<TV_INT> > > regions(TV::m);
    for(int axis=0;axis<TV::m;axis++)Find_Boundary_Regions(regions(axis),mpi_grid.Parallel_Face_Sentinels(axis),false,RANGE<VECTOR<int,1> >(0,0),false);
    ARRAY<MPI_PACKAGE> packages(regions(0).m);ARRAY<ARRAY<T> > buffers(regions(0).m);
    MPI::Datatype T_type=MPI_UTILITIES::Datatype<T>();
    // send and receive into temporary buffers
    int tag=Get_Unique_Tag();
    ARRAY<MPI::Request> requests;
    for(int n=0;n<regions(0).m;n++)if(side_neighbor_ranks(n)!=MPI::PROC_NULL){int axis=n/2;
        packages(n)=mpi_grid.Package_Common_Face_Data(data,axis,regions(axis)(n));
        if(n%2==0) requests.Append(packages(n).Isend(*comm,side_neighbor_ranks(n),tag));
        else requests.Append(packages(n).Irecv(*comm,side_neighbor_ranks(n),tag));}
    // wait
    MPI_UTILITIES::Wait_All(requests);
    MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Assert_Common_Face_Data
//#####################################################################
template<class TV> template<class T_MPI_GRID,class ARRAY<T,FACE_INDEX<TV::m> >> void MPI_GRID<TV>::
Assert_Common_Face_Data(const T_MPI_GRID& mpi_grid,ARRAY<T,FACE_INDEX<TV::m> >& data,const T tolerance) const
{   
    ARRAY<ARRAY<RANGE<TV_INT> > > regions(TV::m);
    for(int axis=0;axis<TV::m;axis++)Find_Boundary_Regions(regions(axis),mpi_grid.Parallel_Face_Sentinels(axis),false,RANGE<VECTOR<int,1> >(0,0),false);
    ARRAY<MPI_PACKAGE> packages(regions(0).m);ARRAY<ARRAY<T> > buffers(regions(0).m);
    MPI::Datatype T_type=MPI_UTILITIES::Datatype<T>();
    // send and receive into temporary buffers
    ARRAY<MPI::Request> requests;
    for(int n=0;n<regions(0).m;n++)if(side_neighbor_ranks(n)!=MPI::PROC_NULL){int axis=n/2;
        packages(n)=mpi_grid.Package_Common_Face_Data(data,axis,regions(axis)(n));
        requests.Append(packages(n).Isend(*comm,side_neighbor_ranks(n),Get_Send_Tag(side_neighbor_directions(n))));
        buffers(n).Resize(packages(n).Size()/sizeof(T));
        requests.Append(comm->Irecv(buffers(n).Get_Array_Pointer(),buffers(n).m,T_type,side_neighbor_ranks(n),Get_Recv_Tag(side_neighbor_directions(n))));}
    // wait
    MPI_UTILITIES::Wait_All(requests);
    // average received data with local data (TODO: find a cleaner general way to do this)
    for(int n=0;n<regions(0).m;n++)if(side_neighbor_ranks(n)!=MPI::PROC_NULL){
        ARRAY<char> pack_buffer(packages(n).Pack_Size(*comm));packages(n).Pack(pack_buffer,*comm);
        ARRAY<T> local_buffer(buffers(n).m);
        int position=0;T_type.Unpack(pack_buffer.Get_Array_Pointer(),pack_buffer.m,local_buffer.Get_Array_Pointer(),local_buffer.m,position,*comm);
        for(int i=0;i<local_buffer.m;i++) assert(abs(buffers(n)(i)-local_buffer(i))<=tolerance);
        position=0;T_type.Pack(local_buffer.Get_Array_Pointer(),local_buffer.m,pack_buffer.Get_Array_Pointer(),pack_buffer.m,position,*comm);
        packages(n).Unpack(pack_buffer,*comm);}
    MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Sync_Common_Face_Weights
//#####################################################################
template<class TV> void MPI_GRID<TV>::
Sync_Common_Face_Weights_From(ARRAY<ARRAY<PAIR<FACE_INDEX<TV::m>,T> >,FACE_INDEX<TV::m> >& weights_to,ARRAY<ARRAY<PAIR<FACE_INDEX<TV::m>,int> >,FACE_INDEX<TV::m> >& weights_from,const int ghost_cells)
{
    int tag=Get_Unique_Tag();

    RANGE<VECTOR<int,2> > range(VECTOR<int,2>(1,1),VECTOR<int,2>(number_of_processes,number_of_processes));
    ARRAY<int,VECTOR<int,2> > all_sizes(range);
    ARRAY<ARRAY<TRIPLE<FACE_INDEX<TV::m>,FACE_INDEX<TV::m>,T> > > transfer_per_proc(number_of_processes); 
    ARRAY<ARRAY<TRIPLE<FACE_INDEX<TV::m>,FACE_INDEX<TV::m>,T> > > receive_per_proc(number_of_processes); 
    for(int axis=0;axis<TV::m;axis++) for(int side=0;side<2;side++){
        int other_rank=side_neighbor_ranks(2*axis+side);if(other_rank<0) continue;
        RANGE<TV_INT> domain=local_grid.Domain_Indices();
        if(side==0) domain.max_corner(axis)=domain.min_corner(axis)+(ghost_cells-1);
        else domain.min_corner(axis)=domain.max_corner(axis)-(ghost_cells-1);
        for(int axis2=0;axis2<TV::m;axis2++){
            RANGE<TV_INT> face_domain=local_grid.Domain_Indices(ghost_cells);
            if(side==0) face_domain.min_corner(axis)+=ghost_cells;
            else face_domain.max_corner(axis)-=ghost_cells;
            if(axis2==axis) face_domain.min_corner(axis)++;
            else face_domain.max_corner(axis2)++;
            for(FACE_ITERATOR<TV> iterator(local_grid,domain,axis2);iterator.Valid();iterator.Next()){FACE_INDEX<TV::m> cell=iterator.Full_Index();
                ARRAY<PAIR<FACE_INDEX<TV::m>,T> >& local_weights=weights_to(cell);
                for(int i=0;i<local_weights.m;i++){
                    if(!face_domain.Lazy_Inside_Half_Open(local_weights(i).x.index)) 
                        transfer_per_proc(other_rank+1).Append(TRIPLE<FACE_INDEX<TV::m>,FACE_INDEX<TV::m>,T>(FACE_INDEX<TV::m>(cell.axis,cell.index+local_to_global_offset),
                            FACE_INDEX<TV::m>(local_weights(i).x.axis,local_weights(i).x.index+local_to_global_offset),local_weights(i).y));}}}}
    for(int i=0;i<number_of_processes;i++) all_sizes(rank+1,i)=transfer_per_proc(i).m;
    for(int i=0;i<number_of_processes;i++) for(int j=0;j<number_of_processes;j++){int output;MPI_UTILITIES::Reduce(all_sizes(i,j),output,MPI::MAX,*comm);all_sizes(i,j)=output;}
    for(int i=0;i<number_of_processes;i++) receive_per_proc(i).Resize(all_sizes(i,rank+1));
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > recv_buffers(number_of_processes);ARRAY<ARRAY<char> > send_buffer(number_of_processes);
    for(int i=0;i<number_of_processes;i++){if(transfer_per_proc(i).m==0) continue;
        int buffer_size=MPI_UTILITIES::Datatype<TV_INT>().Pack_size(transfer_per_proc(i).m*2,*comm)+MPI_UTILITIES::Datatype<int>().Pack_size(transfer_per_proc(i).m*2,*comm)+MPI_UTILITIES::Datatype<T>().Pack_size(transfer_per_proc(i).m,*comm);
        send_buffer(i).Resize(buffer_size);int position=0;
        for(int j=0;j<transfer_per_proc(i).m;j++) MPI_UTILITIES::Pack(transfer_per_proc(i)(j).x.index,transfer_per_proc(i)(j).x.axis,transfer_per_proc(i)(j).y.index,transfer_per_proc(i)(j).y.axis,transfer_per_proc(i)(j).z,send_buffer(i),position,*comm);
        requests.Append(comm->Isend(&send_buffer(i)(0),position,MPI::PACKED,i-1,tag));}
    for(int i=0;i<number_of_processes;i++){if(receive_per_proc(i).m==0) continue;
        MPI::Status status;comm->Probe(i-1,tag,status);
        recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
        requests.Append(comm->Irecv(recv_buffers(i).m?&(recv_buffers(i)(0)):0,recv_buffers(i).m,MPI_PACKED,i-1,tag));}
    MPI_UTILITIES::Wait_All(requests);
    for(int i=0;i<number_of_processes;i++){int position=0;
        for(int j=0;j<receive_per_proc(i).m;j++){
            TRIPLE<FACE_INDEX<TV::m>,FACE_INDEX<TV::m>,T>& data=receive_per_proc(i)(j);
            MPI_UTILITIES::Unpack(data.x.index,data.x.axis,data.y.index,data.y.axis,data.z,recv_buffers(i),position,*comm);
            data.x.index-=local_to_global_offset;data.y.index-=local_to_global_offset;
            RANGE<TV_INT> face_domain=local_grid.Domain_Indices();face_domain.max_corner(data.x.axis)++;
            if(ignore_boundary_faces && face_domain.Lazy_Inside_Half_Open(data.x.index)) continue;
            int index=-1;for(int j=0;j<weights_to(data.x).m;j++) if(weights_to(data.x)(j).x==data.y) index=j;
            assert(data.z>-1e-5);
            if(data.z<0) data.z=0;
            if(index>=0) weights_to(data.x)(index).y=data.z;
            else{
                weights_to(data.x).Append(PAIR<FACE_INDEX<TV::m>,T>(data.y,data.z));
                weights_from(data.y).Append(PAIR<FACE_INDEX<TV::m>,int>(data.x,weights_to(data.x).m));}}}
}
template<class TV> void MPI_GRID<TV>::
Sync_Common_Face_Weights_To(ARRAY<ARRAY<PAIR<FACE_INDEX<TV::m>,T> >,FACE_INDEX<TV::m> >& weights_to,ARRAY<ARRAY<PAIR<FACE_INDEX<TV::m>,int> >,FACE_INDEX<TV::m> >& weights_from,const int ghost_cells)
{
    int tag=Get_Unique_Tag();

    RANGE<VECTOR<int,2> > range(VECTOR<int,2>(1,1),VECTOR<int,2>(number_of_processes,number_of_processes));
    ARRAY<int,VECTOR<int,2> > all_sizes(range);
    ARRAY<ARRAY<TRIPLE<FACE_INDEX<TV::m>,FACE_INDEX<TV::m>,T> > > transfer_per_proc(number_of_processes); 
    ARRAY<ARRAY<TRIPLE<FACE_INDEX<TV::m>,FACE_INDEX<TV::m>,T> > > receive_per_proc(number_of_processes); 
    for(int axis=0;axis<TV::m;axis++) for(int side=0;side<2;side++){
        int other_rank=side_neighbor_ranks(2*axis+side);if(other_rank<0) continue;
        RANGE<TV_INT> domain=local_grid.Domain_Indices();
        if(side==0) domain.max_corner(axis)=domain.min_corner(axis)+(ghost_cells-1);
        else domain.min_corner(axis)=domain.max_corner(axis)-(ghost_cells-1);
        for(int axis2=0;axis2<TV::m;axis2++){
            RANGE<TV_INT> face_domain=local_grid.Domain_Indices(ghost_cells);
            if(side==0) face_domain.min_corner(axis)+=ghost_cells;
            else face_domain.max_corner(axis)-=ghost_cells;
            if(axis2==axis) face_domain.min_corner(axis)++;
            else face_domain.max_corner(axis2)++;
            for(FACE_ITERATOR<TV> iterator(local_grid,domain,axis2);iterator.Valid();iterator.Next()){FACE_INDEX<TV::m> cell=iterator.Full_Index();
                ARRAY<PAIR<FACE_INDEX<TV::m>,int> >& local_weights=weights_from(cell);
                for(int i=0;i<local_weights.m;i++){
                    if(!face_domain.Lazy_Inside_Half_Open(local_weights(i).x.index)) 
                        transfer_per_proc(other_rank+1).Append(TRIPLE<FACE_INDEX<TV::m>,FACE_INDEX<TV::m>,T>(FACE_INDEX<TV::m>(cell.axis,cell.index+local_to_global_offset),
                            FACE_INDEX<TV::m>(local_weights(i).x.axis,local_weights(i).x.index+local_to_global_offset),weights_to(local_weights(i).x)(local_weights(i).y).y));}}}}
    for(int i=0;i<number_of_processes;i++) all_sizes(rank+1,i)=transfer_per_proc(i).m;
    for(int i=0;i<number_of_processes;i++) for(int j=0;j<number_of_processes;j++){int output;MPI_UTILITIES::Reduce(all_sizes(i,j),output,MPI::MAX,*comm);all_sizes(i,j)=output;}
    for(int i=0;i<number_of_processes;i++) receive_per_proc(i).Resize(all_sizes(i,rank+1));
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > recv_buffers(number_of_processes);ARRAY<ARRAY<char> > send_buffer(number_of_processes);
    for(int i=0;i<number_of_processes;i++){if(transfer_per_proc(i).m==0) continue;
        int buffer_size=MPI_UTILITIES::Datatype<TV_INT>().Pack_size(transfer_per_proc(i).m*2,*comm)+MPI_UTILITIES::Datatype<int>().Pack_size(transfer_per_proc(i).m*2,*comm)+MPI_UTILITIES::Datatype<T>().Pack_size(transfer_per_proc(i).m,*comm);
        send_buffer(i).Resize(buffer_size);int position=0;
        for(int j=0;j<transfer_per_proc(i).m;j++) MPI_UTILITIES::Pack(transfer_per_proc(i)(j).x.index,transfer_per_proc(i)(j).x.axis,transfer_per_proc(i)(j).y.index,transfer_per_proc(i)(j).y.axis,transfer_per_proc(i)(j).z,send_buffer(i),position,*comm);
        requests.Append(comm->Isend(&send_buffer(i)(0),position,MPI::PACKED,i-1,tag));}
    for(int i=0;i<number_of_processes;i++){if(receive_per_proc(i).m==0) continue;
        MPI::Status status;comm->Probe(i-1,tag,status);
        recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
        requests.Append(comm->Irecv(recv_buffers(i).m?&(recv_buffers(i)(0)):0,recv_buffers(i).m,MPI_PACKED,i-1,tag));}
    MPI_UTILITIES::Wait_All(requests);
    for(int i=0;i<number_of_processes;i++){int position=0;
        for(int j=0;j<receive_per_proc(i).m;j++){
            TRIPLE<FACE_INDEX<TV::m>,FACE_INDEX<TV::m>,T>& data=receive_per_proc(i)(j);
            MPI_UTILITIES::Unpack(data.x.index,data.x.axis,data.y.index,data.y.axis,data.z,recv_buffers(i),position,*comm);
            data.x.index-=local_to_global_offset;data.y.index-=local_to_global_offset;
            RANGE<TV_INT> face_domain=local_grid.Domain_Indices();face_domain.max_corner(data.x.axis)++;
            if(ignore_boundary_faces && face_domain.Lazy_Inside_Half_Open(data.x.index)) continue;
            int index=-1;for(int j=0;j<weights_to(data.y).m;j++) if(weights_to(data.y)(j).x==data.x) index=j;
            assert(data.z>-1e-5);
            if(data.z<0) data.z=0;
            if(index>=0) weights_to(data.y)(index).y=data.z;
            else{
                weights_to(data.y).Append(PAIR<FACE_INDEX<TV::m>,T>(data.x,data.z));
                weights_from(data.x).Append(PAIR<FACE_INDEX<TV::m>,int>(data.y,weights_to(data.y).m));}}}
}
template<class TV> void MPI_GRID<TV>::
Sync_Common_Cell_Weights_From(ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& weights_to,ARRAY<ARRAY<PAIR<TV_INT,int> >,TV_INT>& weights_from,const int ghost_cells)
{
    int tag=Get_Unique_Tag();
    T delta=1e-4;
    delta=1e-4;

    RANGE<VECTOR<int,2> > range(VECTOR<int,2>(1,1),VECTOR<int,2>(number_of_processes,number_of_processes));
    ARRAY<int,VECTOR<int,2> > all_sizes(range);
    ARRAY<ARRAY<TRIPLE<TV_INT,TV_INT,T> > > transfer_per_proc(number_of_processes); 
    ARRAY<ARRAY<TRIPLE<TV_INT,TV_INT,T> > > receive_per_proc(number_of_processes);
    for(int axis=0;axis<TV::m;axis++) for(int side=0;side<2;side++){
        int other_rank=side_neighbor_ranks(2*axis+side);if(other_rank<0) continue;
        RANGE<TV_INT> domain=local_grid.Domain_Indices(),ghost_domain=local_grid.Domain_Indices(ghost_cells);
        if(side==0) ghost_domain.min_corner(axis)+=ghost_cells;
        else ghost_domain.max_corner(axis)-=ghost_cells;
        if(side==0) domain.max_corner(axis)=domain.min_corner(axis)+(ghost_cells-1);
        else domain.min_corner(axis)=domain.max_corner(axis)-(ghost_cells-1);
        for(CELL_ITERATOR<TV> iterator(local_grid,domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to(cell);
            for(int i=0;i<local_weights.m;i++){
                if(!ghost_domain.Lazy_Inside_Half_Open(local_weights(i).x)) 
                    transfer_per_proc(other_rank+1).Append(TRIPLE<TV_INT,TV_INT,T>(local_to_global_offset+cell,local_to_global_offset+local_weights(i).x,local_weights(i).y));}}}
    for(int i=0;i<number_of_processes;i++) all_sizes(rank+1,i)=transfer_per_proc(i).m;
    for(int i=0;i<number_of_processes;i++) for(int j=0;j<number_of_processes;j++){int output;MPI_UTILITIES::Reduce(all_sizes(i,j),output,MPI::MAX,*comm);all_sizes(i,j)=output;}
    for(int i=0;i<number_of_processes;i++) receive_per_proc(i).Resize(all_sizes(i,rank+1));
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > recv_buffers(number_of_processes);ARRAY<ARRAY<char> > send_buffer(number_of_processes);
    for(int i=0;i<number_of_processes;i++){if(transfer_per_proc(i).m==0) continue;
        int buffer_size=MPI_UTILITIES::Datatype<TV_INT>().Pack_size(transfer_per_proc(i).m*2,*comm)+MPI_UTILITIES::Datatype<T>().Pack_size(transfer_per_proc(i).m,*comm);
        send_buffer(i).Resize(buffer_size);int position=0;
        for(int j=0;j<transfer_per_proc(i).m;j++) MPI_UTILITIES::Pack(transfer_per_proc(i)(j).x,transfer_per_proc(i)(j).y,transfer_per_proc(i)(j).z,send_buffer(i),position,*comm);
        requests.Append(comm->Isend(&send_buffer(i)(0),position,MPI::PACKED,i-1,tag));}
    for(int i=0;i<number_of_processes;i++){if(receive_per_proc(i).m==0) continue;
        MPI::Status status;comm->Probe(i-1,tag,status);
        recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
        requests.Append(comm->Irecv(recv_buffers(i).m?&(recv_buffers(i)(0)):0,recv_buffers(i).m,MPI_PACKED,i-1,tag));}
    MPI_UTILITIES::Wait_All(requests);
    for(int i=0;i<number_of_processes;i++){int position=0;
        for(int j=0;j<receive_per_proc(i).m;j++){
            TRIPLE<TV_INT,TV_INT,T>& data=receive_per_proc(i)(j);
            MPI_UTILITIES::Unpack(data.x,data.y,data.z,recv_buffers(i),position,*comm);
            TV_INT target_cell=data.y-local_to_global_offset;
            int index=-1;for(int k=0;k<weights_to(data.x-local_to_global_offset).m;k++) if(weights_to(data.x-local_to_global_offset)(k).x==target_cell) index=k;
            assert(data.z>-delta);
            if(data.z<0) data.z=0;
            if(index>=0) weights_to(data.x-local_to_global_offset)(index).y=data.z;
            else{
                weights_to(data.x-local_to_global_offset).Append(PAIR<TV_INT,T>(target_cell,data.z));
                weights_from(target_cell).Append(PAIR<TV_INT,int>(data.x-local_to_global_offset,weights_to(data.x-local_to_global_offset).m));}}}
}
template<class TV> void MPI_GRID<TV>::
Sync_Common_Cell_Weights_To(ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& weights_to,ARRAY<ARRAY<PAIR<TV_INT,int> >,TV_INT>& weights_from,const int ghost_cells)
{
    int tag=Get_Unique_Tag();
    RANGE<VECTOR<int,2> > range(VECTOR<int,2>(1,1),VECTOR<int,2>(number_of_processes,number_of_processes));
    ARRAY<int,VECTOR<int,2> > all_sizes(range);
    ARRAY<ARRAY<TRIPLE<TV_INT,TV_INT,T> > > transfer_per_proc(number_of_processes); 
    ARRAY<ARRAY<TRIPLE<TV_INT,TV_INT,T> > > receive_per_proc(number_of_processes);
    for(int axis=0;axis<TV::m;axis++) for(int side=0;side<2;side++){
        int other_rank=side_neighbor_ranks(2*axis+side);if(other_rank<0) continue;
        RANGE<TV_INT> domain=local_grid.Domain_Indices(),ghost_domain=local_grid.Domain_Indices(ghost_cells);
        if(side==0) ghost_domain.min_corner(axis)+=ghost_cells;
        else ghost_domain.max_corner(axis)-=ghost_cells;
        if(side==0) domain.max_corner(axis)=domain.min_corner(axis)+(ghost_cells-1);
        else domain.min_corner(axis)=domain.max_corner(axis)-(ghost_cells-1);
        for(CELL_ITERATOR<TV> iterator(local_grid,domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            ARRAY<PAIR<TV_INT,int> >& local_weights=weights_from(cell);
            for(int i=0;i<local_weights.m;i++)
                if(!ghost_domain.Lazy_Inside_Half_Open(local_weights(i).x)){
                    transfer_per_proc(other_rank+1).Append(TRIPLE<TV_INT,TV_INT,T>(local_to_global_offset+cell,local_to_global_offset+local_weights(i).x,weights_to(local_weights(i).x)(local_weights(i).y).y));}}}
    for(int i=0;i<number_of_processes;i++) all_sizes(rank+1,i)=transfer_per_proc(i).m;
    for(int i=0;i<number_of_processes;i++) for(int j=0;j<number_of_processes;j++){int output;MPI_UTILITIES::Reduce(all_sizes(i,j),output,MPI::MAX,*comm);all_sizes(i,j)=output;}
    for(int i=0;i<number_of_processes;i++) receive_per_proc(i).Resize(all_sizes(i,rank+1));
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > recv_buffers(number_of_processes);ARRAY<ARRAY<char> > send_buffer(number_of_processes);
    for(int i=0;i<number_of_processes;i++){if(transfer_per_proc(i).m==0) continue;
        int buffer_size=MPI_UTILITIES::Datatype<TV_INT>().Pack_size(transfer_per_proc(i).m*2,*comm)+MPI_UTILITIES::Datatype<T>().Pack_size(transfer_per_proc(i).m,*comm);
        send_buffer(i).Resize(buffer_size);int position=0;
        for(int j=0;j<transfer_per_proc(i).m;j++) MPI_UTILITIES::Pack(transfer_per_proc(i)(j).x,transfer_per_proc(i)(j).y,transfer_per_proc(i)(j).z,send_buffer(i),position,*comm);
        requests.Append(comm->Isend(&send_buffer(i)(0),position,MPI::PACKED,i-1,tag));}
    for(int i=0;i<number_of_processes;i++){if(receive_per_proc(i).m==0) continue;
        MPI::Status status;comm->Probe(i-1,tag,status);
        recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
        requests.Append(comm->Irecv(recv_buffers(i).m?&(recv_buffers(i)(0)):0,recv_buffers(i).m,MPI_PACKED,i-1,tag));}
    MPI_UTILITIES::Wait_All(requests);
    for(int i=0;i<number_of_processes;i++){int position=0;
        for(int j=0;j<receive_per_proc(i).m;j++){
            TRIPLE<TV_INT,TV_INT,T>& data=receive_per_proc(i)(j);
            MPI_UTILITIES::Unpack(data.x,data.y,data.z,recv_buffers(i),position,*comm);
            TV_INT target_cell=data.x-local_to_global_offset;
            int index=-1;for(int k=0;k<weights_to(data.y-local_to_global_offset).m;k++) if(weights_to(data.y-local_to_global_offset)(k).x==target_cell) index=k;
            assert(data.z>-1e-5);
            if(data.z<0) data.z=0;
            if(index>=0) weights_to(data.y-local_to_global_offset)(index).y=data.z;
            else{
                weights_to(data.y-local_to_global_offset).Append(PAIR<TV_INT,T>(target_cell,data.z));
                weights_from(target_cell).Append(PAIR<TV_INT,int>(data.y-local_to_global_offset,weights_to(data.y-local_to_global_offset).m));}}}
}
//#####################################################################
// Function Reduce_Add
//#####################################################################
template<class TV> template<class T2> void MPI_GRID<TV>::
Reduce_Add(const T2& input,T2& output) const
{
    MPI_UTILITIES::Reduce(input,output,MPI::SUM,*comm);
}
//#####################################################################
// Function Reduce_Add
//#####################################################################
template<class TV> template<class T2> T2 MPI_GRID<TV>::
Reduce_Add(const T2& local_value) const
{
    T global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::SUM,*comm);
    return global_value;
}
//#####################################################################

#else

//#####################################################################
template<class TV> MPI_GRID<TV>::MPI_GRID(GRID<TV>& local_grid_input,const int number_of_ghost_cells_input,const bool,const TV_INT&,const VECTOR<bool,TV::m>&,MPI::Group* group_input)
    :local_grid(local_grid_input),number_of_ghost_cells(number_of_ghost_cells_input){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> MPI_GRID<TV>::~MPI_GRID(){}
template<class TV> void MPI_GRID<TV>::Initialize_Communicator(const bool manual,MPI::Group* group){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_GRID<TV>::Initialize(VECTOR<VECTOR<bool,2>,TV::m>& domain_walls){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool MPI_GRID<TV>::Neighbor(const int axis,const int axis_side) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_GRID<TV>::Synchronize_Dt(T&) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_GRID<TV>::Synchronize_J_Bounds(int&,int&) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_MPI_GRID,class T2> void MPI_GRID<TV>::Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& data,const int bandwidth,
    const bool include_corners) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_MPI_GRID,class T2> void MPI_GRID<TV>::Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,const GRID<TV>& local_grid,ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& data,
    const int bandwidth,const bool include_corners) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_MPI_GRID,class T2,class T_ARRAYS,class INDEX> void MPI_GRID<TV>::Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,ARRAY_BASE<T2,T_ARRAYS,INDEX>& data,const int bandwidth,
    const bool include_corners) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_MPI_GRID,class T2,class T_ARRAYS,class INDEX> void MPI_GRID<TV>::Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,const GRID<TV>& local_grid,ARRAY_BASE<T2,T_ARRAYS,INDEX>& data,
    const int bandwidth,const bool include_corners) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_MPI_GRID,class T_FACE_ARRAYS> void MPI_GRID<TV>::Exchange_Boundary_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data,const int bandwidth) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_MPI_GRID,class T_FACE_ARRAYS> void MPI_GRID<TV>::Average_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_MPI_GRID,class T_FACE_ARRAYS_BOOL> void MPI_GRID<TV>::Union_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS_BOOL& data) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_MPI_GRID,class T_FACE_ARRAYS> void MPI_GRID<TV>::Copy_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T_MPI_GRID,class T_FACE_ARRAYS> void MPI_GRID<TV>::Assert_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data,const T tolerance) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T2> void MPI_GRID<TV>::Reduce_Add(const T2& input,T2& output) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T2> T2 MPI_GRID<TV>::Reduce_Add(const T2& local_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_GRID<TV>::Sync_Common_Face_Weights_From(ARRAY<ARRAY<PAIR<FACE_INDEX<TV::m>,T> >,FACE_INDEX<TV::m> >& weights_to,ARRAY<ARRAY<PAIR<FACE_INDEX<TV::m>,int> >,FACE_INDEX<TV::m> >& weights_from,const int ghost_cells){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_GRID<TV>::Sync_Common_Face_Weights_To(ARRAY<ARRAY<PAIR<FACE_INDEX<TV::m>,T> >,FACE_INDEX<TV::m> >& weights_to,ARRAY<ARRAY<PAIR<FACE_INDEX<TV::m>,int> >,FACE_INDEX<TV::m> >& weights_from,const int ghost_cells){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_GRID<TV>::Sync_Common_Cell_Weights_From(ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& weights_to,ARRAY<ARRAY<PAIR<TV_INT,int> >,TV_INT>& weights_from,const int ghost_cells) {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_GRID<TV>::Sync_Common_Cell_Weights_To(ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& weights_to,ARRAY<ARRAY<PAIR<TV_INT,int> >,TV_INT>& weights_from,const int ghost_cells) {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################

#endif

//#####################################################################
namespace PhysBAM{
template class MPI_GRID<VECTOR<float,1> >;
template class MPI_GRID<VECTOR<float,2> >;
template class MPI_GRID<VECTOR<float,3> >;
template void MPI_GRID<VECTOR<float,1> >::Assert_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<float,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<float,FACE_INDEX<1> >&,float) const;
template void MPI_GRID<VECTOR<float,1> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<SYMMETRIC_MATRIX<float,1>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<SYMMETRIC_MATRIX<float,1>,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<float,1> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<VECTOR<float,1>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<VECTOR<float,1>,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<float,1> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<VECTOR<float,3>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<VECTOR<float,3>,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<float,1> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<float,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<float,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<float,1> >::Copy_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<float,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<float,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,SYMMETRIC_MATRIX<float,1> >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAYS_ND_BASE<SYMMETRIC_MATRIX<float,1>,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,VECTOR<float,1> >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAYS_ND_BASE<VECTOR<float,1>,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,VECTOR<float,3> >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAYS_ND_BASE<VECTOR<float,3>,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,GRID<VECTOR<float,1> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,float>(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAYS_ND_BASE<float,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<SYMMETRIC_MATRIX<float,1>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<SYMMETRIC_MATRIX<float,1>,FACE_INDEX<1> >&,int) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<VECTOR<float,1>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<VECTOR<float,1>,FACE_INDEX<1> >&,int) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<VECTOR<float,3>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<VECTOR<float,3>,FACE_INDEX<1> >&,int) const;
template void MPI_GRID<VECTOR<float,1> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<float,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<float,FACE_INDEX<1> >&,int) const;
template void MPI_GRID<VECTOR<float,1> >::Reduce_Add<ARRAY<float,int> >(ARRAY<float,int> const&,ARRAY<float,int>&) const;
template void MPI_GRID<VECTOR<float,1> >::Union_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,1> >,ARRAY<bool,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,ARRAY<bool,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<float,2> >::Assert_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<float,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<float,FACE_INDEX<2> >&,float) const;
template void MPI_GRID<VECTOR<float,2> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<SYMMETRIC_MATRIX<float,2>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<SYMMETRIC_MATRIX<float,2>,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<float,2> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<VECTOR<float,2>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<VECTOR<float,2>,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<float,2> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<VECTOR<float,4>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<VECTOR<float,4>,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<float,2> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<float,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<float,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<float,2> >::Copy_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<float,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<float,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,SYMMETRIC_MATRIX<float,2> >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAYS_ND_BASE<SYMMETRIC_MATRIX<float,2>,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,VECTOR<float,2> >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAYS_ND_BASE<VECTOR<float,2>,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,VECTOR<float,4> >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAYS_ND_BASE<VECTOR<float,4>,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,GRID<VECTOR<float,2> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,float>(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAYS_ND_BASE<float,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<SYMMETRIC_MATRIX<float,2>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<SYMMETRIC_MATRIX<float,2>,FACE_INDEX<2> >&,int) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<VECTOR<float,2>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<VECTOR<float,2>,FACE_INDEX<2> >&,int) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<VECTOR<float,4>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<VECTOR<float,4>,FACE_INDEX<2> >&,int) const;
template void MPI_GRID<VECTOR<float,2> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<float,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<float,FACE_INDEX<2> >&,int) const;
template void MPI_GRID<VECTOR<float,2> >::Reduce_Add<ARRAY<float,int> >(ARRAY<float,int> const&,ARRAY<float,int>&) const;
template void MPI_GRID<VECTOR<float,2> >::Union_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,2> >,ARRAY<bool,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,ARRAY<bool,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<float,3> >::Assert_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<float,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<float,FACE_INDEX<3> >&,float) const;
template void MPI_GRID<VECTOR<float,3> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<SYMMETRIC_MATRIX<float,3>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<SYMMETRIC_MATRIX<float,3>,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<float,3> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<VECTOR<float,3>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<VECTOR<float,3>,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<float,3> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<VECTOR<float,5>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<VECTOR<float,5>,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<float,3> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<float,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<float,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<float,3> >::Copy_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<float,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<float,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,SYMMETRIC_MATRIX<float,3> >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAYS_ND_BASE<SYMMETRIC_MATRIX<float,3>,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,VECTOR<float,3> >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAYS_ND_BASE<VECTOR<float,3>,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,VECTOR<float,5> >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAYS_ND_BASE<VECTOR<float,5>,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,GRID<VECTOR<float,3> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,float>(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAYS_ND_BASE<float,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<SYMMETRIC_MATRIX<float,3>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<SYMMETRIC_MATRIX<float,3>,FACE_INDEX<3> >&,int) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<VECTOR<float,3>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<VECTOR<float,3>,FACE_INDEX<3> >&,int) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<VECTOR<float,5>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<VECTOR<float,5>,FACE_INDEX<3> >&,int) const;
template void MPI_GRID<VECTOR<float,3> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<float,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<float,FACE_INDEX<3> >&,int) const;
template void MPI_GRID<VECTOR<float,3> >::Reduce_Add<ARRAY<float,int> >(ARRAY<float,int> const&,ARRAY<float,int>&) const;
template void MPI_GRID<VECTOR<float,3> >::Union_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<float,3> >,ARRAY<bool,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,ARRAY<bool,FACE_INDEX<3> >&) const;
template class MPI_GRID<VECTOR<double,1> >;
template class MPI_GRID<VECTOR<double,2> >;
template class MPI_GRID<VECTOR<double,3> >;
template void MPI_GRID<VECTOR<double,1> >::Assert_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<double,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<double,FACE_INDEX<1> >&,double) const;
template void MPI_GRID<VECTOR<double,1> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<SYMMETRIC_MATRIX<double,1>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<SYMMETRIC_MATRIX<double,1>,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<double,1> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<VECTOR<double,1>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<VECTOR<double,1>,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<double,1> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<VECTOR<double,3>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<VECTOR<double,3>,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<double,1> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<double,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<double,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<double,1> >::Copy_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<double,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<double,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,SYMMETRIC_MATRIX<double,1> >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAYS_ND_BASE<SYMMETRIC_MATRIX<double,1>,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,VECTOR<double,1> >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAYS_ND_BASE<VECTOR<double,1>,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,VECTOR<double,3> >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAYS_ND_BASE<VECTOR<double,3>,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,GRID<VECTOR<double,1> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,double>(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAYS_ND_BASE<double,VECTOR<int,1> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<SYMMETRIC_MATRIX<double,1>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<SYMMETRIC_MATRIX<double,1>,FACE_INDEX<1> >&,int) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<VECTOR<double,1>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<VECTOR<double,1>,FACE_INDEX<1> >&,int) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<VECTOR<double,3>,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<VECTOR<double,3>,FACE_INDEX<1> >&,int) const;
template void MPI_GRID<VECTOR<double,1> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<double,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<double,FACE_INDEX<1> >&,int) const;
template void MPI_GRID<VECTOR<double,1> >::Reduce_Add<ARRAY<double,int> >(ARRAY<double,int> const&,ARRAY<double,int>&) const;
template void MPI_GRID<VECTOR<double,1> >::Union_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,1> >,ARRAY<bool,FACE_INDEX<1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,ARRAY<bool,FACE_INDEX<1> >&) const;
template void MPI_GRID<VECTOR<double,2> >::Assert_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<double,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<2> >&,double) const;
template void MPI_GRID<VECTOR<double,2> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<SYMMETRIC_MATRIX<double,2>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<SYMMETRIC_MATRIX<double,2>,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<double,2> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<VECTOR<double,2>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<VECTOR<double,2>,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<double,2> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<VECTOR<double,4>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<VECTOR<double,4>,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<double,2> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<double,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<double,2> >::Copy_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<double,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,SYMMETRIC_MATRIX<double,2> >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAYS_ND_BASE<SYMMETRIC_MATRIX<double,2>,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,VECTOR<double,2> >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAYS_ND_BASE<VECTOR<double,2>,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,VECTOR<double,4> >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAYS_ND_BASE<VECTOR<double,4>,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,GRID<VECTOR<double,2> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,double>(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAYS_ND_BASE<double,VECTOR<int,2> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<SYMMETRIC_MATRIX<double,2>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<SYMMETRIC_MATRIX<double,2>,FACE_INDEX<2> >&,int) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<VECTOR<double,2>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<VECTOR<double,2>,FACE_INDEX<2> >&,int) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<VECTOR<double,4>,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<VECTOR<double,4>,FACE_INDEX<2> >&,int) const;
template void MPI_GRID<VECTOR<double,2> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<double,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<2> >&,int) const;
template void MPI_GRID<VECTOR<double,2> >::Reduce_Add<ARRAY<double,int> >(ARRAY<double,int> const&,ARRAY<double,int>&) const;
template void MPI_GRID<VECTOR<double,2> >::Union_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,2> >,ARRAY<bool,FACE_INDEX<2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,ARRAY<bool,FACE_INDEX<2> >&) const;
template void MPI_GRID<VECTOR<double,3> >::Assert_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<double,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<double,FACE_INDEX<3> >&,double) const;
template void MPI_GRID<VECTOR<double,3> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<SYMMETRIC_MATRIX<double,3>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<SYMMETRIC_MATRIX<double,3>,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<double,3> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<VECTOR<double,3>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<VECTOR<double,3>,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<double,3> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<VECTOR<double,5>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<VECTOR<double,5>,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<double,3> >::Average_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<double,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<double,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<double,3> >::Copy_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<double,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<double,FACE_INDEX<3> >&) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,SYMMETRIC_MATRIX<double,3> >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAYS_ND_BASE<SYMMETRIC_MATRIX<double,3>,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,VECTOR<double,3> >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAYS_ND_BASE<VECTOR<double,3>,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,VECTOR<double,5> >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAYS_ND_BASE<VECTOR<double,5>,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,bool>(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,GRID<VECTOR<double,3> > const&,ARRAYS_ND_BASE<bool,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Cell_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,double>(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAYS_ND_BASE<double,VECTOR<int,3> >&,int,bool) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<SYMMETRIC_MATRIX<double,3>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<SYMMETRIC_MATRIX<double,3>,FACE_INDEX<3> >&,int) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<VECTOR<double,3>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<VECTOR<double,3>,FACE_INDEX<3> >&,int) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<VECTOR<double,5>,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<VECTOR<double,5>,FACE_INDEX<3> >&,int) const;
template void MPI_GRID<VECTOR<double,3> >::Exchange_Boundary_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<double,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<double,FACE_INDEX<3> >&,int) const;
template void MPI_GRID<VECTOR<double,3> >::Reduce_Add<ARRAY<double,int> >(ARRAY<double,int> const&,ARRAY<double,int>&) const;
template void MPI_GRID<VECTOR<double,3> >::Union_Common_Face_Data<MPI_UNIFORM_GRID<VECTOR<double,3> >,ARRAY<bool,FACE_INDEX<3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,ARRAY<bool,FACE_INDEX<3> >&) const;
}
