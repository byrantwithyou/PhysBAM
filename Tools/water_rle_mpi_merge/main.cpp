//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_SIMPLE_ITERATOR.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#include <PhysBAM_Tools/Particles/Delete_Particles.h>
using namespace PhysBAM;
//#####################################################################
// Class LOCAL_GRID
//#####################################################################
template<class T_GRID>
class LOCAL_GRID
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;
    typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;typedef typename T_GRID::HORIZONTAL_GRID T_HORIZONTAL_GRID;

    T_GRID grid;
    MPI_RLE_GRID<T_GRID> mpi_grid;
    const T_GRID& global_grid;
    const typename GRID<TV>& global_uniform_grid;
    TV_HORIZONTAL_INT offset;
    ARRAY<bool> neighbor_overlaps;

    LOCAL_GRID(const T_GRID& global_grid_input,const typename GRID<TV>& global_uniform_grid_input)
        :mpi_grid(grid,true),global_grid(global_grid_input),global_uniform_grid(global_uniform_grid_input)
    {}

    void Initialize()
    {// find offset
    TV_HORIZONTAL X=grid.horizontal_grid.Domain().Minimum_Corner();
    TV_HORIZONTAL_INT local_offset,global_offset;
    local_offset=grid.horizontal_grid.Closest_Node(X);
    global_offset=global_uniform_grid.Get_Horizontal_Grid().Closest_Node(X);
    offset=global_offset-local_offset;
    // pretend there are neighbors on all sides for use in Find_Boundary_Regions
    mpi_grid.side_neighbor_ranks.Resize(T_HORIZONTAL_GRID::number_of_neighbors_per_cell);mpi_grid.side_neighbor_ranks.Fill(1);
    mpi_grid.all_neighbor_ranks.Resize(T_HORIZONTAL_GRID::number_of_one_ring_neighbors_per_cell);mpi_grid.all_neighbor_ranks.Fill(1);
    ARRAY<RANGE<TV_HORIZONTAL_INT> > regions;mpi_grid.Find_Boundary_Regions(regions,RANGE<TV_HORIZONTAL_INT>::Zero_Box(),true,RANGE<VECTOR<int,1> >(-5,-5),true);
    RANGE<TV_HORIZONTAL_INT> global_region=global_uniform_grid.Get_Horizontal_Grid().Domain_Indices();
    neighbor_overlaps.Resize(regions.m);
    for(int n=0;n<regions.m;n++)if(regions(n).Lazy_Intersection(global_region-offset)) neighbor_overlaps(n)=true;}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,grid);Initialize();}

    RANGE<TV_HORIZONTAL_INT> Interior_Region(const RANGE<TV_HORIZONTAL_INT>& sentinels) const
    {return grid.horizontal_grid.Get_MAC_Grid().Domain_Indices()+sentinels;}

    template<class T_ITERATOR,class T2> void Put(const ARRAY<T2>& local_data,const RANGE<TV_HORIZONTAL_INT>& region,ARRAY<T2>& global_data) const
    {RLE_GRID_SIMPLE_ITERATOR<T_GRID,T_ITERATOR> local(grid,region),global(global_grid,region+offset);
    for(;local;local++,global++)global_data(global.index)=local_data(local.index);}

    template<class T_ITERATOR,class T2> void Put(const ARRAY<T2>& local_data,ARRAY<T2>& global_data) const
    {Put<T_ITERATOR>(local_data,Interior_Region(T_ITERATOR::Sentinels()),global_data);
    ARRAY<RANGE<TV_HORIZONTAL_INT> > regions;mpi_grid.Find_Boundary_Regions(regions,T_ITERATOR::Sentinels(),true,RANGE<VECTOR<int,1> >(-grid.number_of_ghost_cells,-1),true);
    for(int n=0;n<regions.m;n++)if(!neighbor_overlaps(n)) Put<T_ITERATOR>(local_data,regions(n),global_data);}

    template<class T2> struct Put_Component{template<class T_FACE> static void Apply(const LOCAL_GRID& local_grid,const ARRAY<T2>& local_data,ARRAY<T2>& global_data)
    {local_grid.Put<T_FACE>(local_data,global_data);}};

    template<class T_ITERATOR,class T2> void Get(ARRAY<T2>& local_data,const RANGE<TV_HORIZONTAL_INT>& region,const ARRAY<T2>& global_data) const
    {RLE_GRID_SIMPLE_ITERATOR<T_GRID,T_ITERATOR> local(grid,region),global(global_grid,region+offset);
    for(;local;local++,global++)local_data(local.index)=global_data(global.index);}

    template<class T_ITERATOR,class T2> void Get(ARRAY<T2>& local_data,const ARRAY<T2>& global_data) const
    {Get<T_ITERATOR>(local_data,Interior_Region(T_ITERATOR::Sentinels()).Thickened(grid.number_of_ghost_cells),global_data);}

    template<class T2> struct Get_Component{template<class T_FACE> static void Apply(const LOCAL_GRID& local_grid,ARRAY<T2>& local_data,const ARRAY<T2>& global_data)
    {local_grid.Get<T_FACE>(local_data,global_data);}};

    template<class T_ITERATOR,class T2> T Maximum_Error(const ARRAY<T2>& local_data,const ARRAY<T2>& global_data,const int bandwidth,int& index) const
    {RANGE<TV_HORIZONTAL_INT> region=Interior_Region(T_ITERATOR::Sentinels()).Thickened(bandwidth);
    T max_error=0;RLE_GRID_SIMPLE_ITERATOR<T_GRID,T_ITERATOR> local(grid,region),global(global_grid,region+offset);
    for(;local;local++,global++){T error=(T)fabs(global_data(global.index)-local_data(local.index));
        if(max_error<error){max_error=error;index=local.index;}}
    return max_error;}
};
//#####################################################################
// Class RLE_MERGER
//#####################################################################
template<class T_GRID,class RW>
class RLE_MERGER
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename GRID<TV> T_UNIFORM_GRID;typedef typename T_GRID::FACE_Y_ITERATOR FACE_Y_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::BLOCK_ITERATOR BLOCK_ITERATOR;
    typedef typename T_GRID::RUN T_RUN;typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;typedef typename T_GRID::ARRAYS_HORIZONTAL T_ARRAYS_HORIZONTAL;
    typedef typename T_GRID::HORIZONTAL_GRID T_HORIZONTAL_GRID;typedef typename T_HORIZONTAL_GRID::CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;
    typedef typename T_HORIZONTAL_GRID::NODE_ITERATOR HORIZONTAL_NODE_ITERATOR;

    int number_of_processes;
    TV_HORIZONTAL_INT processes_per_dimension;
    const std::string input_directory,output_directory;
    T_GRID grid;
    ARRAY<LOCAL_GRID<T_GRID>*> local_grids;
    T_UNIFORM_GRID global_uniform_grid;
    bool merge;

    bool process_levelset;
    bool process_object_levelset;
    bool process_debug_data;
    bool process_particles;
    bool process_removed_particles;
    bool process_velocities;

    RLE_MERGER(const std::string& input_directory_input,const std::string& output_directory_input,PARSE_ARGS& parse_args)
        :input_directory(input_directory_input),output_directory(output_directory_input)
    {
        FILE_UTILITIES::Read_From_File<RW>(input_directory+"/1/global_uniform_grid",global_uniform_grid);
        process_levelset=!parse_args.Get_Option_Value("-skip_levelset");
        process_object_levelset=!parse_args.Get_Option_Value("-skip_object_levelset");
        process_debug_data=!parse_args.Get_Option_Value("-skip_debug_data");
        process_particles=!parse_args.Get_Option_Value("-skip_particles");
        process_removed_particles=!parse_args.Get_Option_Value("-skip_removed_particles");
        process_velocities=!parse_args.Get_Option_Value("-skip_velocities");
        if(parse_args.Get_Option_Value("-fast")) process_object_levelset=process_debug_data=process_particles=process_removed_particles=process_velocities=false;
        if(parse_args.Get_Option_Value("-removed_particles")) process_removed_particles=true;
        processes_per_dimension=Parse_Split_Option(parse_args,TV_HORIZONTAL_INT());
        if(processes_per_dimension==TV_HORIZONTAL_INT()){merge=true;Get_Number_Of_Processes(parse_args);}
        else merge=false;
    }

    void Get_Number_Of_Processes(PARSE_ARGS& parse_args)
    {number_of_processes=parse_args.Get_Integer_Value("-np");
    if(number_of_processes<0){LOG::cerr<<"Invalid np<0"<<std::endl;exit(1);}
    else if(!number_of_processes){ // autodetect number of processes
        FILE_UTILITIES::Find_First_Nonexistent_Directory_In_Sequence(STRING_UTILITIES::string_sprintf("%s/%%d",input_directory),1,&number_of_processes);--number_of_processes;
        LOG::cout<<"Autodetected "<<number_of_processes<<" processes"<<std::endl;}}

    VECTOR<int,1> Parse_Split_Option(PARSE_ARGS& parse_args,const VECTOR<int,1>&)
    {if(parse_args.Is_Value_Set("-split")) PHYSBAM_FATAL_ERROR();
    return VECTOR<int,1>();}

    VECTOR<int,2> Parse_Split_Option(PARSE_ARGS& parse_args,const VECTOR<int,2>&)
    {VECTOR<int,2> number_of_process_per_dimension(parse_args.Get_Vector_2D_Value("-split"));
    number_of_processes=number_of_process_per_dimension.x*number_of_process_per_dimension.y;
    return number_of_process_per_dimension;}

    bool Need_Merge(const std::string &filename)
    {std::string output_filename=output_directory+"/"+filename;
    if(!FILE_UTILITIES::File_Exists(output_filename)) return true;
    // check file times
    std::string prefix=input_directory+"/"+filename+".";
    for(int i=0;i<number_of_processes;i++)if(FILE_UTILITIES::Compare_File_Times(input_directory+STRING_UTILITIES::string_sprintf("/%d/",i)+filename,output_filename)>0) return true;
    return false;}

    bool Source_Files_Exist(const int frame) const
    {std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
        for(int i=0;i<number_of_processes;i++)if(!FILE_UTILITIES::File_Exists(input_directory+"/"+STRING_UTILITIES::string_sprintf("%d",i)+"/rle_grid."+f)) return false;
    return true;}

    static bool Grid_Parameters_Match(const T_GRID& grid1,const T_GRID& grid2)
    {return grid1.number_of_ghost_cells==grid2.number_of_ghost_cells && grid1.minimum_vertical_space==grid2.minimum_vertical_space && grid1.minimum_long_run_length==grid2.minimum_long_run_length
        && grid1.long_run_cells==grid2.long_run_cells && grid1.long_run_faces_horizontal==grid2.long_run_faces_horizontal && grid1.negative_bandwidth==grid2.negative_bandwidth
        && grid1.positive_bandwidth==grid2.positive_bandwidth && grid1.jmin==grid2.jmin && grid1.jmax==grid2.jmax;}

    template<class T2> struct Horizontal_Maximum_Error{template<class T_FACE> static void Apply(const std::string& prefix,const LOCAL_GRID<T_GRID>& local_grid,
        const ARRAY<T2>& local_data,const ARRAY<T2>& global_data,const int bandwidth,const T threshold=1e-7)
    {int index=0;
    T max_error=local_grid.template Maximum_Error<T_FACE>(local_data,global_data,bandwidth,index);
    if(max_error>threshold){LOG::cout<<prefix<<", face "<<(T_FACE::Direction().x?"x":"z")<<" = "<<max_error<<" ("<<index<<")"<<std::endl;}}};

    template<class T2> void Vertical_Maximum_Error(const std::string& prefix,const LOCAL_GRID<T_GRID>& local_grid,const ARRAY<T2>& local_data,const ARRAY<T2>& global_data,
        const int bandwidth,const T threshold=1e-7)
    {int index=0;
    T max_error=local_grid.template Maximum_Error<FACE_Y_ITERATOR>(local_data,global_data,bandwidth,index);
    if(max_error>threshold){LOG::cout<<prefix<<", face y = "<<max_error<<" ("<<index<<")"<<std::endl;}}

    void Merge_All_Frames(const int first_frame,const int last_frame)
    {for(int frame=first_frame;frame<=last_frame;frame++){
        if(!Source_Files_Exist(frame)) break;
        Merge(frame);}}

//#####################################################################
    void Merge(const int frame);
    void Merge_Grids(const std::string& filename);
    template<class T2> void Merge_Cell_Data(const std::string& filename,const int verify_bandwidth);
    template<class T2> void Merge_Face_Data(const std::string& filename,const int verify_bandwidth);
    template<class T_PARTICLES> void Merge_Particles(const std::string& filename);
//#####################################################################
};
//#####################################################################
// Function Merge
//#####################################################################
template<class T_GRID,class RW> void RLE_MERGER<T_GRID,RW>::
Merge(const int frame)
{
    LOG::SCOPE scope("FRAME","Frame %d",frame);
    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);
    Merge_Grids("rle_grid"+f);
    if(process_levelset) Merge_Cell_Data<T>("rle_levelset"+f,grid.number_of_ghost_cells);
    if(process_object_levelset) Merge_Cell_Data<T>("rle_object_levelset"+f,grid.number_of_ghost_cells);
    if(process_debug_data){
        Merge_Cell_Data<T>("rle_pressure"+f,1);
        Merge_Cell_Data<int>("rle_colors"+f,0); // TODO: consider changing this back to 1
        Merge_Cell_Data<bool>("rle_psi_D"+f,1);
        Merge_Face_Data<bool>("rle_psi_N"+f,grid.number_of_ghost_cells);}
    if(process_velocities) Merge_Face_Data<T>("rle_face_velocities"+f,grid.number_of_ghost_cells);
    if(process_particles){
        Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >("negative_particles"+f);
        Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >("positive_particles"+f);}
    if(process_removed_particles){
        Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >("removed_negative_particles"+f);
        Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >("removed_positive_particles"+f);}
    if(merge) FILE_UTILITIES::Write_To_Text_File(output_directory+"/last_frame",frame,"\n");
    local_grids.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Merge_Grids
//#####################################################################
template<class T_GRID,class RW> void RLE_MERGER<T_GRID,RW>::
Merge_Grids(const std::string& filename)
{
    local_grids.Resize(number_of_processes);
    if(merge){
        for(int p=0;p<number_of_processes;p++){
            local_grids(p)=new LOCAL_GRID<T_GRID>(grid,global_uniform_grid);
            FILE_UTILITIES::Read_From_File<RW>(input_directory+STRING_UTILITIES::string_sprintf("/%d/",p)+filename,*local_grids(p));}

        // verify that grid parameters match
        for(int p=2;p<=number_of_processes;p++)if(!Grid_Parameters_Match(local_grids(1)->grid,local_grids(p)->grid)){
            LOG::cerr<<"Grid parameter mismatch between processes 1 and "<<p<<std::endl;assert(false);exit(1);}

        // merge grids
        T_GRID::Transfer_Noncolumn_Data(local_grids(1)->grid,grid);
        grid.uniform_grid=global_uniform_grid;grid.horizontal_grid=grid.uniform_grid.Get_Horizontal_Grid();

        // sanity check
        int global_column_count=grid.horizontal_grid.Domain_Indices().Size(),total_local_column_count=0;
        for(int p=0;p<number_of_processes;p++) total_local_column_count+=local_grids(p)->grid.horizontal_grid.Domain_Indices().Size();
        if(global_column_count!=total_local_column_count){
            LOG::cout<<"Column mismatch:"<<std::endl;
            LOG::cout<<"  global horizontal grid = "<<grid.horizontal_grid<<", columns = "<<global_column_count<<std::endl;
            for(int p=0;p<number_of_processes;p++)
                LOG::cout<<"  local horizontal grid "<<p<<" = "<<local_grids(p)->grid.horizontal_grid<<", columns = "<<local_grids(p)->grid.horizontal_grid.Domain_Indices().Size()<<std::endl;
            exit(1);}

        // copying columns
        grid.columns.Resize(grid.horizontal_grid.Get_MAC_Grid().Domain_Indices(grid.number_of_ghost_cells+1),false,false);
        for(int p=0;p<number_of_processes;p++)T_GRID::ARRAYS_HORIZONTAL::template REBIND<ARRAY<T_RUN> >::TYPE::Shifted_Put(local_grids(p)->grid.columns,grid.columns,local_grids(p)->offset);
        grid.Topology_Changed();grid.Compute_Auxiliary_Information();
        // verifying column equality
        for(int p=0;p<number_of_processes;p++){
            TV_HORIZONTAL_INT offset=local_grids(p)->offset;
            typename T_ARRAYS_HORIZONTAL::template REBIND<ARRAY<T_RUN> >::TYPE& local_columns=local_grids(p)->grid.columns;
            for(HORIZONTAL_CELL_ITERATOR iterator(local_grids(p)->grid.horizontal_grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
                TV_HORIZONTAL_INT cell=iterator.Cell_Index();
                if(grid.columns(cell+offset)!=local_columns(cell)){
                    LOG::cerr<<"Interior column equality check failed"<<std::endl;
                    LOG::cerr<<"  process "<<p<<", column = "<<cell<<std::endl;
                    assert(false);exit(1);}}
            for(HORIZONTAL_CELL_ITERATOR iterator(local_grids(p)->grid.horizontal_grid,grid.number_of_ghost_cells+1);iterator.Valid();iterator.Next()){
                TV_HORIZONTAL_INT cell=iterator.Cell_Index();
                if(grid.columns(cell+offset)!=local_columns(cell)){
                    LOG::cerr<<"Far ghost column equality check failed"<<std::endl;
                    LOG::cerr<<"  process "<<p<<", column = "<<cell<<std::endl;}}}
        FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,grid);}
    else{
        FILE_UTILITIES::Read_From_File<RW>(output_directory+"/"+filename,grid);
        MPI_GRID<T_HORIZONTAL_GRID> mpi_grid(grid.horizontal_grid,true);
        mpi_grid.Split_Grid(processes_per_dimension);
        for(HORIZONTAL_NODE_ITERATOR iterator(mpi_grid.process_grid);iterator.Valid();iterator.Next()){
            TV_HORIZONTAL_INT coordinates=iterator.Node_Index();
            int p=mpi_grid.process_ranks(coordinates);
            local_grids(p)=new LOCAL_GRID<T_GRID>(grid,grid.uniform_grid);
            local_grids(p)->mpi_grid.Update_Local_Grid(mpi_grid.Restrict_Grid(coordinates));
            T_GRID::Transfer_Noncolumn_Data(grid,local_grids(p)->grid);
            local_grids(p)->Initialize();
            // copy columns
            local_grids(p)->grid.columns.Resize(local_grids(p)->grid.horizontal_grid.Get_MAC_Grid().Domain_Indices(grid.number_of_ghost_cells+1),false,false);
            T_GRID::ARRAYS_HORIZONTAL::template REBIND<ARRAY<T_RUN> >::TYPE::Shifted_Get(local_grids(p)->grid.columns,grid.columns,-local_grids(p)->offset);
            local_grids(p)->grid.Topology_Changed();local_grids(p)->grid.Compute_Auxiliary_Information();
            FILE_UTILITIES::Write_To_File<RW>(input_directory+STRING_UTILITIES::string_sprintf("/%d/",p)+filename,local_grids(p)->grid);}}
}
//#####################################################################
// Function Merge_Cell_Data
//#####################################################################
template<class T_GRID,class RW> template<class T2> void RLE_MERGER<T_GRID,RW>::
Merge_Cell_Data(const std::string& filename,const int verify_bandwidth)
{
    if(merge){
        // read
        ARRAY<ARRAY<T2> > local_data(number_of_processes);
        for(int p=0;p<number_of_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",p)+filename;
            if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return;}
            FILE_UTILITIES::Read_From_File<RW>(name,local_data(p));}
        // merge
        ARRAY<T2> global_data(grid.number_of_cells);
        for(int p=0;p<number_of_processes;p++)local_grids(p)->template Put<CELL_ITERATOR>(local_data(p),global_data);
        // verify
        for(int p=0;p<number_of_processes;p++){int index=0;
            T max_error=local_grids(p)->template Maximum_Error<CELL_ITERATOR>(local_data(p),global_data,verify_bandwidth,index);
            if(max_error>0){LOG::cout<<filename<<": max error on process "<<p<<" = "<<max_error<<" ("<<index<<" = "<<local_data(p)(index)<<")"<<std::endl;}}
        // write
        FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);}
    else{
        std::string name=output_directory+"/"+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping split"<<std::endl;return;}
        ARRAY<T2> global_data;FILE_UTILITIES::Read_From_File<RW>(name,global_data);
        for(int p=0;p<number_of_processes;p++){
            ARRAY<T2> local_data(local_grids(p)->grid.number_of_cells);
            local_grids(p)->template Get<CELL_ITERATOR>(local_data,global_data);
            FILE_UTILITIES::Write_To_File<RW>(input_directory+STRING_UTILITIES::string_sprintf("/%d/",p)+filename,local_data);}}
}
//#####################################################################
// Function Merge_Face_Data
//#####################################################################
template<class T_GRID,class RW> template<class T2> void RLE_MERGER<T_GRID,RW>::
Merge_Face_Data(const std::string& filename,const int verify_bandwidth)
{
    if(merge){
        // read
        ARRAY<ARRAY<T2> > local_data(number_of_processes);
        for(int p=0;p<number_of_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",p)+filename;
            if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return;}
            FILE_UTILITIES::Read_From_File<RW>(name,local_data(p));}
        // merge
        ARRAY<T2> global_data(grid.number_of_faces);
        for(int p=0;p<number_of_processes;p++)T_GRID::template Face_Loop<typename LOCAL_GRID<T_GRID>::template Put_Component<T2> >(*local_grids(p),local_data(p),global_data);
        // verify
        for(int p=0;p<number_of_processes;p++){
            std::string prefix=filename+": max error on process "+STRING_UTILITIES::string_sprintf("%d",p);
            T_GRID::template Horizontal_Face_Loop<Horizontal_Maximum_Error<T2> >(prefix,*local_grids(p),local_data(p),global_data,verify_bandwidth);
            Vertical_Maximum_Error<T2>(prefix,*local_grids(p),local_data(p),global_data,verify_bandwidth);}
        // write
        FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);}
    else{
        std::string name=output_directory+"/"+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping split"<<std::endl;return;}
        ARRAY<T2> global_data;FILE_UTILITIES::Read_From_File<RW>(name,global_data);
        for(int p=0;p<number_of_processes;p++){
            ARRAY<T2> local_data(local_grids(p)->grid.number_of_cells);
            T_GRID::template Face_Loop<typename LOCAL_GRID<T_GRID>::template Get_Component<T2> >(*local_grids(p),local_data,global_data);
            FILE_UTILITIES::Write_To_File<RW>(input_directory+STRING_UTILITIES::string_sprintf("/%d/",p)+filename,local_data);}}
}
//#####################################################################
// Function Merge_Particles
//#####################################################################
template<class T_GRID,class RW> template<class T_PARTICLES> void RLE_MERGER<T_GRID,RW>::
Merge_Particles(const std::string& filename)
{
    if(merge){
        // read
        int process_without_particles=0;
        bool particles_in_long_cells=false;
        ARRAY<ARRAY<T_PARTICLES*> > local_data(number_of_processes);
        for(int p=0;p<number_of_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",p)+filename;
            if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return;}
            FILE_UTILITIES::Read_From_File<RW>(name,local_data(p));
            if(local_data(p).m!=local_grids(p)->grid.number_of_blocks && local_data(p).m!=local_grids(p)->grid.number_of_blocks+1){
                if(!local_data(p).m) process_without_particles=p;
                else{LOG::cerr<<"filename"<<": process "<<p<<" has incorrectly sized particle array."<<std::endl;exit(1);}}
            if(local_data(p).m==local_grids(p)->grid.number_of_blocks+1) particles_in_long_cells=true;}
        // check for zero size
        if(process_without_particles){
            for(int p=0;p<number_of_processes;p++)if(local_data(p).m){
                LOG::cerr<<filename<<": process "<<p<<" has particles but process "<<process_without_particles<<" does not."<<std::endl;exit(1);}
            FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,ARRAY<T_PARTICLES*>());return;}
        // merge
        ARRAY<T_PARTICLES*> global_data(grid.number_of_blocks+particles_in_long_cells);
        for(int p=0;p<number_of_processes;p++){
            RANGE<TV_HORIZONTAL_INT> region=local_grids(p)->Interior_Region(BLOCK_ITERATOR::Sentinels()),interior_region=region.Thickened(-1);
            // copy interior particles
            {BLOCK_ITERATOR local(local_grids(p)->grid,interior_region),global(grid,interior_region+local_grids(p)->offset);
            for(;local;local++,global++)exchange(global_data(global.Block()),local_data(p)(local.Block()));}
            // merge boundary particles
            {BLOCK_ITERATOR local(local_grids(p)->grid,region),global(grid,region+local_grids(p)->offset);
            for(;local;local++,global++){
                T_PARTICLES *&from_particles=local_data(p)(local.Block()),*&to_particles=global_data(global.Block());
                if(!from_particles) continue;
                if(!to_particles){exchange(to_particles,from_particles);continue;}
                to_particles->Append(*from_particles);}
            // copy particles in long cells
            if(local_data(p).m==local_grids(p)->grid.number_of_blocks+1 && local_data(p).Last()){
                if(global_data.Last()) global_data.Last()->Take(*local_data(p).Last());
                else{global_data.Last()=local_data(p).Last();local_data(p).Last()=0;}}}}
        // write
        FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
        // free memory
        for(int p=0;p<number_of_processes;p++)local_data(p).Delete_Pointers_And_Clean_Memory();
        global_data.Delete_Pointers_And_Clean_Memory();}
    else{
        std::string name=output_directory+"/"+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping split"<<std::endl;return;}
        ARRAY<T_PARTICLES*> global_data;FILE_UTILITIES::Read_From_File<RW>(name,global_data);
        if(!global_data.m) return;
        bool particles_in_long_cells=global_data.m==grid.number_of_blocks+1;
        for(int p=0;p<number_of_processes;p++){
            ARRAY<T_PARTICLES*> local_data(local_grids(p)->grid.number_of_blocks+particles_in_long_cells);
            RANGE<TV_HORIZONTAL_INT> region=local_grids(p)->Interior_Region(BLOCK_ITERATOR::Sentinels()),interior_region=region.Thickened(-1);
            // copy interior particles
            {BLOCK_ITERATOR local(local_grids(p)->grid,interior_region),global(grid,interior_region+local_grids(p)->offset);
            for(;local;local++,global++)exchange(global_data(global.Block()),local_data(local.Block()));}
            // split boundary particles
            {BLOCK_ITERATOR local(local_grids(p)->grid,region),global(grid,region+local_grids(p)->offset);
            for(;local;local++,global++){
                T_PARTICLES *&to_particles=local_data(local.Block()),*&from_particles=global_data(global.Block());
                if(!from_particles) continue;
                if(!to_particles) to_particles=from_particles->Clone();
                else to_particles->Initialize(*from_particles);
                Delete_Particles_Outside_Grid(*to_particles,local_grids(p)->grid);}}
            // split particles in long cells
            if(particles_in_long_cells){
                T_PARTICLES *&to_particles=local_data(local_data.m),*&from_particles=global_data(global_data.m);
                if(!from_particles) continue;
                if(!to_particles) to_particles=from_particles->Clone();
                else to_particles->Initialize(*from_particles);
                Delete_Particles_Outside_Grid(*to_particles,local_grids(p)->grid);}
            // write
            FILE_UTILITIES::Write_To_File<RW>(input_directory+STRING_UTILITIES::string_sprintf("/%d/",p)+filename,local_data);
            local_data.Delete_Pointers_And_Clean_Memory();}
        global_data.Delete_Pointers_And_Clean_Memory();}
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-f","force");
    parse_args.Add_Option_Argument("-2d","assume 2d");
    parse_args.Add_Option_Argument("-3d","assume 3d");
    parse_args.Add_Integer_Argument("-start_frame",0,"start frame number");
    parse_args.Add_Integer_Argument("-np",0,"number of mpi processes (use 0 to autodetect)");
    parse_args.Add_String_Argument("-o","","output directory");
    parse_args.Add_Vector_2D_Argument("-split",VECTOR<double,2>(0,0),"split into given number of processes per dimension");
    parse_args.Add_Option_Argument("-skip_levelset","skip_levelset");
    parse_args.Add_Option_Argument("-skip_object_levelset","skip_object_levelset");
    parse_args.Add_Option_Argument("-skip_debug_data","skip_debug_data");
    parse_args.Add_Option_Argument("-skip_particles","skip_particles");
    parse_args.Add_Option_Argument("-skip_removed_particles","skip_removed_particles");
    parse_args.Add_Option_Argument("-skip_velocities","skip_velocities");
    parse_args.Add_Option_Argument("-fast","skip everything except for levelset");
    parse_args.Add_Option_Argument("-removed_particles","process removed particles");
    parse_args.Set_Extra_Arguments(-1,"<input_directory>");
    parse_args.Parse(argc,argv);

    std::string input_directory;
    if(!parse_args.Num_Extra_Args()) input_directory=".";
    else if(parse_args.Num_Extra_Args()==1) input_directory=parse_args.Extra_Arg(1);
    else parse_args.Print_Usage(true);
    std::string output_directory=input_directory;
    if(parse_args.Is_Value_Set("-o")) output_directory=parse_args.Get_String_Value("-o");

    int first_frame,last_frame;
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/1/first_frame",first_frame);
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/1/last_frame",last_frame);
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/first_frame",first_frame);

    if(parse_args.Is_Value_Set("-start_frame")) first_frame=parse_args.Get_Integer_Value("-start_frame");

    bool assume_2d=parse_args.Get_Option_Value("-2d"),assume_3d=parse_args.Get_Option_Value("-3d");
    if(!assume_2d && !assume_3d){ // guess dimension
        std::string grid_file=input_directory+"/1/global_uniform_grid";
        if(!FILE_UTILITIES::File_Exists(grid_file)) grid_file=input_directory+"/rle_grid.0";
        VECTOR<int,3> mnmn;FILE_UTILITIES::Read_From_File<float>(input_directory+"/1/global_uniform_grid",mnmn);
        if(mnmn.z>10 && mnmn.z<20000) assume_3d=true;else assume_2d=true;}

    if(assume_2d){
        RLE_MERGER<RLE_GRID_2D<float>,float> merger(input_directory,output_directory,parse_args);
        merger.Merge_All_Frames(first_frame,last_frame);}
    else{
        RLE_MERGER<RLE_GRID_3D<float>,float> merger(input_directory,output_directory,parse_args);
        merger.Merge_All_Frames(first_frame,last_frame);}
    LOG::cout<<std::endl;
}
//#####################################################################
