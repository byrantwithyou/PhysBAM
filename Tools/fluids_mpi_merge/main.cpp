//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/LOCAL_GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Class MERGER
//#####################################################################
template<class T,class T_GRID,class RW>
class MERGER
{
public:
    typedef typename T_GRID::VECTOR_T TV;typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS;typedef UNIFORM_GRID_ITERATOR_NODE<TV> NODE_ITERATOR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<VECTOR<T,T_GRID::dimension+2> >::TYPE T_ARRAYS_DIMENSION_SCALAR;
    

    const int number_of_processes;
    int number_of_fluid_processes;
    int fluid_proc_offset;
    const std::string input_directory,output_directory;
    T_GRID grid;
    ARRAY<LOCAL_GRID<T_GRID>*> local_grids;

    ARRAY<int> needs_init,needs_destroy;
    SOLID_BODY_COLLECTION<TV>* solid_body_collection;

    bool merge_levelset;
    bool merge_object_levelset;
    bool merge_debug_data;
    bool merge_particles;
    bool merge_removed_particles;
    bool merge_removed_particle_times;
    bool merge_velocities;
    bool merge_temperatures;
    bool merge_densities;
    bool merge_compressible;
    bool merge_solid_fluid;
    bool print_max_errors;

    MERGER(const int number_of_processes_input,const std::string& input_directory_input,const std::string& output_directory_input,const bool print_max_errors,const bool merge_solid_fluid_input)
        :number_of_processes(number_of_processes_input),number_of_fluid_processes(number_of_processes_input),fluid_proc_offset(0),input_directory(input_directory_input),
        output_directory(output_directory_input),merge_levelset(true),merge_object_levelset(true),merge_debug_data(true),merge_particles(true),merge_removed_particles(true),
        merge_removed_particle_times(false),merge_velocities(true),merge_temperatures(true),merge_densities(true),merge_compressible(false),merge_solid_fluid(merge_solid_fluid_input),
        print_max_errors(print_max_errors)
    {
        if(merge_solid_fluid){
            number_of_fluid_processes=number_of_processes-1;
            fluid_proc_offset=number_of_processes-number_of_fluid_processes;
            FILE_UTILITIES::Read_From_File<RW>(input_directory+"/2/common/global_grid",grid);}
        else FILE_UTILITIES::Read_From_File<RW>(input_directory+"/1/common/global_grid",grid);
        grid=grid.Get_MAC_Grid();
    }

    ~MERGER()
    {delete solid_body_collection;}

    bool Need_Merge(const std::string &filename)
    {std::string output_filename=output_directory+"/"+filename;
    if(!FILE_UTILITIES::File_Exists(output_filename)) return true;
    // check file times
    std::string prefix=input_directory+"/"+filename+".";
    for(int i=0;i<number_of_processes;i++)if(FILE_UTILITIES::Compare_File_Times(input_directory+STRING_UTILITIES::string_sprintf("/%d/",i)+filename,output_filename)>0) return true;
    return false;}

    bool Source_Files_Exist(const int frame) const
    {for(int i=0;i<number_of_processes;i++)if(!FILE_UTILITIES::File_Exists(input_directory+"/"+STRING_UTILITIES::string_sprintf("%d/%d",i,frame)+"/time")) return false;
    return true;}

    void Merge_All_Frames(const int first_frame,const int last_frame)
    {
        if(merge_solid_fluid) solid_body_collection=new SOLID_BODY_COLLECTION<TV>(0);
        for(int frame=first_frame;frame<=last_frame;frame++){
            //if(!Source_Files_Exist(frame)){LOG::cout<<"missing source files for frame "<<frame<<std::endl;break;}
        Merge(frame);}
    }

//#####################################################################
    void Merge(const int frame);
    template<class T_ARRAYS> bool Merge_Cell_Data(const std::string& filename,const int verify_bandwidth,const bool scale=false);
    template<class T_ARRAYS> bool Merge_Levelset(const std::string& filename,const int verify_bandwidth);
    template<class T_FACE_ARRAYS_2> bool Merge_Face_Data(const std::string& filename,const int verify_bandwidth);
    template<class T_LIST_2> bool Merge_Lists(const std::string& filename);
    template<class T_PARTICLES> bool Merge_Particles(const std::string& filename);
    template<class T_ARRAYS> void Scale_Cell_Data(T_ARRAYS& array){PHYSBAM_NOT_IMPLEMENTED();}
    void Scale_Cell_Data(T_ARRAYS_SCALAR& array);
//#####################################################################
};
//#####################################################################
// Function Merge
//#####################################################################
template<class T,class T_GRID,class RW> void MERGER<T,T_GRID,RW>::
Merge(const int frame)
{
    LOG::SCOPE scope("FRAME","Frame %d",frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);

    local_grids.Resize(number_of_fluid_processes);
    for(int p=0;p<number_of_fluid_processes;p++){
        local_grids(p)=new LOCAL_GRID<T_GRID>(grid);
        FILE_UTILITIES::Read_From_File<RW>(input_directory+STRING_UTILITIES::string_sprintf("/%d/common/grid",p+fluid_proc_offset),*local_grids(p));}
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid",grid);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    if(merge_compressible){
        Merge_Cell_Data<T_ARRAYS_DIMENSION_SCALAR>(f+"euler_U",0);
        Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"density",0);
        Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"pressure",0);
        Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"temperature",0);
        Merge_Cell_Data<ARRAY<TV,TV_INT> >(f+"centered_velocities",0);
        if(FILE_UTILITIES::File_Exists(input_directory+"/2/"+f+"soot")) Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"soot",0);
        if(FILE_UTILITIES::File_Exists(input_directory+"/2/"+f+"soot_fuel")) Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"soot_fuel",0);
        if(merge_debug_data){
            Merge_Cell_Data<ARRAY<bool,TV_INT> >(f+"psi_D",1);
            Merge_Face_Data<T_FACE_ARRAYS_BOOL>(f+"psi_N",1);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"density_gradient",0,true);
            Merge_Cell_Data<ARRAY<bool,TV_INT> >(f+"euler_psi",1);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"energy",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"enthalpy",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"entropy",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"internal_energy",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"machnumber",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"speedofsound",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"velocity_plus_c",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"velocity_minus_c",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"pressure_gradient",0,true);}
        T time;FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/1/%d/time",input_directory.c_str(),frame),time);
        FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/%d/time",output_directory.c_str(),frame),time);
        std::string levelset_file=input_directory+"/1"+f+"levelset";
        if(FILE_UTILITIES::File_Exists(levelset_file)){
            Merge_Levelset<T_ARRAYS_SCALAR>(f+"levelset",3);
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/levelset_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Levelset<T_ARRAYS_SCALAR>(filename,3)) break;}}}
    else{
        if(merge_levelset){
            Merge_Levelset<T_ARRAYS_SCALAR>(f+"levelset",3);
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/levelset_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Levelset<T_ARRAYS_SCALAR>(filename,3)) break;}}
        if(merge_object_levelset) Merge_Levelset<T_ARRAYS_SCALAR>(f+"object_levelset",3);
        if(merge_debug_data){
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"pressure",1);
            Merge_Cell_Data<T_ARRAYS_INT>(f+"colors",0); // TODO: consider changing this back to 1
            Merge_Cell_Data<ARRAY<bool,TV_INT> >(f+"psi_D",1);
            Merge_Face_Data<T_FACE_ARRAYS_BOOL>(f+"psi_N",1);
            Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >(f+"negative_particles");}
        if(merge_velocities) Merge_Face_Data<T_FACE_ARRAYS>(f+"mac_velocities",3);
        if(merge_temperatures) Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"temperature",3);
        Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"reaction_speed",3);
        if(merge_densities) Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"density",3);
        if(merge_particles){
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/negative_particles_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >(filename)) break;}
            Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >(f+"positive_particles");
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/positive_particles_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >(filename)) break;}
            Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(f+"removed_negative_particles");}
        if(merge_removed_particles){
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/removed_negative_particles_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(filename)) break;}
            Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(f+"removed_positive_particles");
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/removed_positive_particles_%d.%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(filename)) break;}}
        if(merge_removed_particle_times) Merge_Lists<ARRAY<PAIR<int,T> > >(f+"removed_particle_times");}
    if(merge_solid_fluid){
        solid_body_collection->rigid_body_collection.Read(STREAM_TYPE(RW()),STRING_UTILITIES::string_sprintf("%s/1/",input_directory.c_str()),frame,&needs_init);
        solid_body_collection->rigid_body_collection.rigid_geometry_collection.structure_list.Fill_Needs_Write();
        solid_body_collection->rigid_body_collection.Write(STREAM_TYPE(RW()),output_directory,frame);
        bool include_static_frame=frame==0;
        solid_body_collection->deformable_body_collection.Read(STREAM_TYPE(RW()),STRING_UTILITIES::string_sprintf("%s/1/",input_directory.c_str()),frame,-1,include_static_frame,false);
        solid_body_collection->deformable_body_collection.Write(STREAM_TYPE(RW()),output_directory+"/",frame,-1,include_static_frame,false);}
    local_grids.Delete_Pointers_And_Clean_Memory();
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
// Function Merge_Cell_Data
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_ARRAYS> bool MERGER<T,T_GRID,RW>::
Merge_Cell_Data(const std::string& filename,const int verify_bandwidth,const bool scale)
{
    // read
    ARRAY<T_ARRAYS> local_data(number_of_fluid_processes);
    for(int p=0;p<number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        FILE_UTILITIES::Read_From_File<RW>(name,local_data(p));}
    // merge
    T_ARRAYS global_data(grid.Cell_Indices(3));
    for(int p=0;p<number_of_fluid_processes;p++)local_grids(p)->Put(local_data(p),global_data);
    // verify
    for(int p=0;p<number_of_fluid_processes;p++){TV_INT index;
        T max_error=local_grids(p)->Maximum_Error(local_data(p),global_data,verify_bandwidth,index);
        if(max_error>0 && print_max_errors){LOG::cout<<filename<<": max error on process "<<p<<" = "<<max_error<<" ("<<index<<" = "<<local_data(p)(index)<<")"<<std::endl;}}
    if(scale) Scale_Cell_Data(global_data);
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    return true;
}
//#####################################################################
// Function Scale_Cell_Data
//#####################################################################
template<class T,class T_GRID,class RW> void MERGER<T,T_GRID,RW>::
Scale_Cell_Data(T_ARRAYS_SCALAR& array)
{
    T max_val=array.Max();
    if(max_val)
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next())
            array(iterator.Cell_Index())=exp(-array(iterator.Cell_Index())/max_val);
}
//#####################################################################
// Function Merge_Levelset
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_ARRAYS> bool MERGER<T,T_GRID,RW>::
Merge_Levelset(const std::string& filename,const int verify_bandwidth)
{
    // read
    ARRAY<T_ARRAYS> local_data(number_of_fluid_processes);
    for(int p=0;p<number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        T_GRID skip;FILE_UTILITIES::Read_From_File<RW>(name,skip,local_data(p));}
    // merge
    T_ARRAYS global_data(grid.Cell_Indices(3));
    for(int p=0;p<number_of_fluid_processes;p++)local_grids(p)->Put(local_data(p),global_data);
    // verify
    for(int p=0;p<number_of_fluid_processes;p++){TV_INT index;
        T max_error=local_grids(p)->Maximum_Error(local_data(p),global_data,verify_bandwidth,index);
        if(max_error>0 && print_max_errors){LOG::cout<<filename<<": max error on process "<<p<<" = "<<max_error<<" ("<<index<<" = "<<local_data(p)(index)<<")"<<std::endl;}}
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,grid,global_data);
    return true;
}
//#####################################################################
// Function Merge_Face_Data
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_FACE_ARRAYS_2> bool MERGER<T,T_GRID,RW>::
Merge_Face_Data(const std::string& filename,const int verify_bandwidth)
{
    // read
    ARRAY<T_FACE_ARRAYS_2*> local_data(number_of_fluid_processes);
    for(int p=0;p<number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        local_data(p)=new T_FACE_ARRAYS_2;
        FILE_UTILITIES::Read_From_File<RW>(name,*local_data(p));}
    // merge
    T_FACE_ARRAYS_2 global_data(grid,3);
    for(int p=0;p<number_of_fluid_processes;p++)local_grids(p)->Put_Faces(*local_data(p),global_data);
    // verify
    for(int p=0;p<number_of_fluid_processes;p++){
        if(print_max_errors){
            std::string prefix=filename+": max error on process "+STRING_UTILITIES::string_sprintf("%d",p);
            local_grids(p)->Maximum_Error(prefix,*local_data(p),global_data,verify_bandwidth);}}
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    local_data.Delete_Pointers_And_Clean_Memory();
    return true;
}
//#####################################################################
// Function Merge_Lists
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_LIST_2> bool MERGER<T,T_GRID,RW>::
Merge_Lists(const std::string& filename)
{
    // read
    ARRAY<T_LIST_2*> local_data(number_of_fluid_processes);
    for(int p=0;p<number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        local_data(p)=new T_LIST_2;
        FILE_UTILITIES::Read_From_File<RW>(name,*local_data(p));}
    // merge
    T_LIST_2 global_data;
    for(int p=0;p<number_of_fluid_processes;p++) global_data.Append_Elements(*local_data(p));
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    local_data.Delete_Pointers_And_Clean_Memory();
    return true;
}
//#####################################################################
// Function Merge_Particles
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_PARTICLES> bool MERGER<T,T_GRID,RW>::
Merge_Particles(const std::string& filename)
{
    // read
    int process_without_particles=0;
    ARRAY<typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE> local_data(number_of_fluid_processes);
    for(int p=0;p<number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        FILE_UTILITIES::Read_From_File<RW>(name,local_data(p));}
    // check for zero size
    if(process_without_particles){
        for(int p=0;p<number_of_fluid_processes;p++)if(local_data(p).Size().x){
            LOG::cerr<<filename<<": process "<<p<<" has particles but process "<<process_without_particles<<" does not."<<std::endl;exit(1);}
        FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE());return true;}
    // merge
    typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE global_data(grid.Node_Indices());
    for(int p=0;p<number_of_fluid_processes;p++){
        RANGE<TV_INT> region=local_grids(p)->Interior_Region(RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector())),interior_region=region.Thickened(-1);
        // copy interior particles
        {NODE_ITERATOR local(local_grids(p)->grid,interior_region),global(grid,interior_region+local_grids(p)->offset);
        for(;local.Valid();local.Next(),global.Next())exchange(global_data(global.Node_Index()),local_data(p)(local.Node_Index()));}
        // merge boundary particles
        {NODE_ITERATOR local(local_grids(p)->grid,region),global(grid,region+local_grids(p)->offset);
        for(;local.Valid();local.Next(),global.Next()){
            T_PARTICLES *&from_particles=local_data(p)(local.Node_Index()),*&to_particles=global_data(global.Node_Index());
            if(!from_particles) continue;
            if(!to_particles){exchange(to_particles,from_particles);continue;}
            to_particles->Append(*from_particles);}}}
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    // free memory
    for(int p=0;p<number_of_fluid_processes;p++)local_data(p).Delete_Pointers_And_Clean_Memory();
    global_data.Delete_Pointers_And_Clean_Memory();
    return true;
}
//#####################################################################
// Function Do_Merge
//#####################################################################
template<class T> void
Do_Merge(PARSE_ARGS& parse_args)
{

    int first_frame=0,last_frame=0,number_of_processes=0;
    bool opt_minimal=false,print_maxerrors=false,is_solid_fluid=false;
    bool is_1d=false,is_2d=false,is_3d=false;
    bool opt_merge_levelset=false,opt_merge_object_levelset=false,opt_merge_debug_data=false;
    bool opt_merge_particles=false,opt_merge_removed_particles=false,opt_merge_removed_particle_times=false;
    bool opt_merge_velocities=false,opt_merge_compressible=false;
    std::string input_directory,output_directory;
    parse_args.Add("-start_frame",&first_frame,"frame","start frame number");
    parse_args.Add("-last_frame",&last_frame,"frame","last frame number");
    parse_args.Add("-np",&number_of_processes,"number","number of mpi processes (use 0 to autodetect)");
    parse_args.Add("-o",&output_directory,"dir","output directory");
    parse_args.Add_Not("-skip_levelset",&opt_merge_levelset,"skip_levelset");
    parse_args.Add_Not("-skip_object_levelset",&opt_merge_object_levelset,"skip_object_levelset");
    parse_args.Add_Not("-skip_debug_data",&opt_merge_debug_data,"skip_debug_data");
    parse_args.Add_Not("-skip_particles",&opt_merge_particles,"skip_particles");
    parse_args.Add_Not("-skip_removed_particles",&opt_merge_removed_particles,"skip_removed_particles");
    parse_args.Add("-removed_particle_times",&opt_merge_removed_particle_times,"removed_particle_times");
    parse_args.Add_Not("-skip_velocities",&opt_merge_velocities,"skip_velocities");
    parse_args.Add("-minimal",&opt_minimal,"skip everything but the levelset");
    parse_args.Add("-print_maxerrors",&print_maxerrors,"print max errors");
    parse_args.Add("-compressible",&opt_merge_compressible,"input data is compressible output");
    parse_args.Add("-solid_fluid",&is_solid_fluid,"input data is solid fluid data");
    parse_args.Add("-1d",&is_1d,"input data is 1-D");
    parse_args.Add("-2d",&is_2d,"input data is 2-D");
    parse_args.Add("-3d",&is_3d,"input data is 3-D");
    parse_args.Extra(&input_directory,"input_directory","input_directory");
    parse_args.Parse();

    FILE_UTILITIES::Read_From_Text_File(input_directory+"/1/common/first_frame",first_frame);
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/1/common/last_frame",last_frame);
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/first_frame",first_frame);

    if(number_of_processes<0){LOG::cerr<<"Invalid np<0"<<std::endl;exit(1);}
    else if(!number_of_processes){ // autodetect number of processes
        FILE_UTILITIES::Find_First_Nonexistent_Directory_In_Sequence(STRING_UTILITIES::string_sprintf("%s/%%d",input_directory.c_str()),1,&number_of_processes);--number_of_processes;
        LOG::cout<<"Autodetected "<<number_of_processes<<" processes"<<std::endl;}

    if(!(is_1d || is_2d || is_3d)) {
        VECTOR<int,3> mnmn;
        if(!is_solid_fluid) FILE_UTILITIES::Read_From_File<T>(input_directory+"/1/common/grid",mnmn);
        else FILE_UTILITIES::Read_From_File<T>(input_directory+"/2/common/grid",mnmn);
        if(mnmn.z>10 && mnmn.z<20000) is_3d=true;
        else if(mnmn.y>10 && mnmn.y<20000) is_2d=true;
        else is_1d=true;}

    if(is_1d){
        MERGER<T,GRID<VECTOR<T,1> >,T> merger(number_of_processes,input_directory,output_directory,print_maxerrors,is_solid_fluid);
        merger.merge_levelset=opt_merge_levelset;
        merger.merge_object_levelset=opt_merge_object_levelset;
        merger.merge_debug_data=opt_merge_debug_data;
        merger.merge_particles=opt_merge_particles;
        merger.merge_removed_particles=opt_merge_removed_particles;
        merger.merge_removed_particle_times=opt_merge_removed_particle_times;
        merger.merge_velocities=opt_merge_velocities;
        merger.merge_compressible=opt_merge_compressible;
        if(opt_minimal){
            merger.merge_object_levelset=false;merger.merge_debug_data=false;merger.merge_particles=false;merger.merge_removed_particles=false;merger.merge_removed_particle_times=false;merger.merge_velocities=false;}
        merger.Merge_All_Frames(first_frame,last_frame);}
    else if(is_2d){
        MERGER<T,GRID<VECTOR<T,2> >,T> merger(number_of_processes,input_directory,output_directory,print_maxerrors,is_solid_fluid);
        merger.merge_levelset=opt_merge_levelset;
        merger.merge_object_levelset=opt_merge_object_levelset;
        merger.merge_debug_data=opt_merge_debug_data;
        merger.merge_particles=opt_merge_particles;
        merger.merge_removed_particles=opt_merge_removed_particles;
        merger.merge_removed_particle_times=opt_merge_removed_particle_times;
        merger.merge_velocities=opt_merge_velocities;
        merger.merge_compressible=opt_merge_compressible;
        if(opt_minimal){
            merger.merge_object_levelset=false;merger.merge_debug_data=false;merger.merge_particles=false;merger.merge_removed_particles=false;merger.merge_removed_particle_times=false;merger.merge_velocities=false;}
        merger.Merge_All_Frames(first_frame,last_frame);}
    else if(is_3d){
        MERGER<T,GRID<VECTOR<T,3> >,T> merger(number_of_processes,input_directory,output_directory,print_maxerrors,is_solid_fluid);
        merger.merge_levelset=opt_merge_levelset;
        merger.merge_object_levelset=opt_merge_object_levelset;
        merger.merge_debug_data=opt_merge_debug_data;
        merger.merge_particles=opt_merge_particles;
        merger.merge_removed_particles=opt_merge_removed_particles;
        merger.merge_removed_particle_times=opt_merge_removed_particle_times;
        merger.merge_velocities=opt_merge_velocities;
        merger.merge_compressible=opt_merge_compressible;
        if(opt_minimal){
            merger.merge_object_levelset=false;merger.merge_debug_data=false;merger.merge_particles=false;merger.merge_removed_particles=false;merger.merge_removed_particle_times=false;merger.merge_velocities=false;}
        merger.Merge_All_Frames(first_frame,last_frame);}
    LOG::cout<<std::endl;
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char* argv[])
{
    bool opt_double=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-double",&opt_double,"input data is in doubles");
    parse_args.Parse(true);

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(opt_double) Do_Merge<double>(parse_args); 
    else Do_Merge<float>(parse_args);
#else
    Do_Merge<float>(parse_args);
#endif
}
//#####################################################################
