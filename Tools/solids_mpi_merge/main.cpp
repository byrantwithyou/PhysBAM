//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Log/LOG.h>
#include <Tools/Parallel_Computation/LOCAL_GRID.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Class MERGER
//#####################################################################
template<class T,class TV,class RW>
class MERGER
{
public:
    const int number_of_processes;
    const std::string input_directory,output_directory;
    RIGID_BODY_COLLECTION<TV> rigid_body_collection;
    ARRAY<RIGID_BODY_COLLECTION<TV>*> local_rigid_collections;

    ARRAY<int> needs_init,needs_destroy;
    ARRAY<int> colors;

    bool merge_rigids;
    bool print_max_errors;

    MERGER(const int number_of_processes_input,const std::string& input_directory_input,const std::string& output_directory_input,const bool print_max_errors)
        :number_of_processes(number_of_processes_input),input_directory(input_directory_input),output_directory(output_directory_input),rigid_body_collection(0),merge_rigids(true),print_max_errors(print_max_errors)
    {
    }

    ~MERGER()
    {}

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
        for(int frame=first_frame;frame<=last_frame;frame++){Merge(frame);}
    }

//#####################################################################
    void Merge(const int frame);
    bool Merge_Rigid_Data(const int frame);
    template<class T_LIST_2> bool Merge_Lists(const std::string& filename);
//#####################################################################
};
//#####################################################################
// Function Merge
//#####################################################################
template<class T,class TV,class RW> void MERGER<T,TV,RW>::
Merge(const int frame)
{
    LOG::SCOPE scope("FRAME","Frame %d",frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    if(merge_rigids) Merge_Rigid_Data(frame);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+f+"/colors",colors);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
// Function Merge_Rigid_Data
//#####################################################################
template<class T,class TV,class RW> bool MERGER<T,TV,RW>::
Merge_Rigid_Data(const int frame)
{    
    //int current_index=0;
    local_rigid_collections.Resize(number_of_processes);
    for(int p=0;p<number_of_processes;p++){
        local_rigid_collections(p)=new RIGID_BODY_COLLECTION<TV>(0);
        local_rigid_collections(p)->Read(STREAM_TYPE(RW()),STRING_UTILITIES::string_sprintf("%s/%d/",input_directory.c_str(),p),frame,&needs_init);
        local_rigid_collections(p)->structure_list.Fill_Needs_Write();
        rigid_body_collection.rigid_body_particles.Add_Elements(local_rigid_collections(p)->rigid_body_particles.Size());
        colors.Resize(local_rigid_collections(p)->rigid_body_particles.Size());
        for(int i=0;i<local_rigid_collections(p)->rigid_body_particles.Size();i++){
            if(!local_rigid_collections(p)->rigid_body_particles.rigid_body(i)) continue;
            if(!local_rigid_collections(p)->Is_Active(i)) continue;
            colors(i)=p;
            rigid_body_collection.rigid_body_particles.Copy_Element(local_rigid_collections(p)->rigid_body_particles,i,i);
            rigid_body_collection.Add_Rigid_Body_And_Geometry(&local_rigid_collections(p)->Rigid_Body(i));}}
    rigid_body_collection.Write(STREAM_TYPE(RW()),output_directory,frame);
    return true;
}
//#####################################################################
// Function Merge_Lists
//#####################################################################
template<class T,class TV,class RW> template<class T_LIST_2> bool MERGER<T,TV,RW>::
Merge_Lists(const std::string& filename)
{
    // read
    ARRAY<T_LIST_2*> local_data(number_of_processes);
    for(int p=0;p<number_of_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        local_data(p)=new T_LIST_2;
        FILE_UTILITIES::Read_From_File<RW>(name,*local_data(p));}
    // merge
    T_LIST_2 global_data;
    for(int p=0;p<number_of_processes;p++) global_data.Append_Elements(*local_data(p));
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    local_data.Delete_Pointers_And_Clean_Memory();
    return true;
}
//#####################################################################
// Function Do_Merge
//#####################################################################
template<class T> void
Do_Merge(PARSE_ARGS& parse_args)
{
    std::string output_directory,input_directory;
    bool is_1d=false,is_2d=false,is_3d=false,print_maxerrors=false;
    int first_frame=0,last_frame=0,number_of_processes=0;
    parse_args.Add("-start_frame",&first_frame,"frame","start frame number");
    parse_args.Add("-last_frame",&last_frame,"frame","last frame number");
    parse_args.Add("-np",&number_of_processes,"number","number of mpi processes (use 0 to autodetect)");
    parse_args.Add("-o",&output_directory,"file","output directory");
    parse_args.Add("-print_maxerrors",&print_maxerrors,"print max errors");
    parse_args.Add("-1d",&is_1d,"input data is 1-D");
    parse_args.Add("-2d",&is_2d,"input data is 2-D");
    parse_args.Add("-3d",&is_3d,"input data is 3-D");
    parse_args.Extra(&input_directory,"input_directory","input_directory");
    parse_args.Parse();

    if(!output_directory.size()) output_directory=input_directory+"/merged";

    FILE_UTILITIES::Read_From_Text_File(input_directory+"/1/common/first_frame",first_frame);
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/1/common/last_frame",last_frame);
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/first_frame",first_frame);

    if(!is_1d && !is_2d && !is_3d){LOG::cout<<"Assuming 3d"<<std::endl;is_3d=true;}
    if(number_of_processes<0){LOG::cerr<<"Invalid np<0"<<std::endl;exit(1);}
    else if(!number_of_processes){ // autodetect number of processes
        FILE_UTILITIES::Find_First_Nonexistent_Directory_In_Sequence(STRING_UTILITIES::string_sprintf("%s/%%d",input_directory.c_str()),1,&number_of_processes);--number_of_processes;
        LOG::cout<<"Autodetected "<<number_of_processes<<" processes"<<std::endl;}

    if(is_1d){
        MERGER<T,VECTOR<T,1>,T> merger(number_of_processes,input_directory,output_directory,print_maxerrors);
        merger.Merge_All_Frames(first_frame,last_frame);}
    else if(is_2d){
        MERGER<T,VECTOR<T,2>,T> merger(number_of_processes,input_directory,output_directory,print_maxerrors);
        merger.Merge_All_Frames(first_frame,last_frame);}
    else if(is_3d){
        MERGER<T,VECTOR<T,3>,T> merger(number_of_processes,input_directory,output_directory,print_maxerrors);
        merger.Merge_All_Frames(first_frame,last_frame);}
    LOG::cout<<std::endl;
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char* argv[])
{
    bool use_doubles=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-double",&use_doubles,"input data is in doubles");
    parse_args.Parse(true);

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(use_doubles) Do_Merge<double>(parse_args); 
    else Do_Merge<float>(parse_args);
#else
    Do_Merge<float>(parse_args);
#endif
}
//#####################################################################
