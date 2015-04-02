//#####################################################################
// Copyright 2004, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG_ENTRY.h>
#include <Tools/Utilities/TIMER.h>
namespace PhysBAM{
namespace LOG{
//#####################################################################
bool LOG_ENTRY::start_on_separate_line=false;
bool LOG_ENTRY::log_file_start_on_separate_line=false;
bool LOG_ENTRY::needs_indent=true;
bool LOG_ENTRY::log_file_needs_indent=true;
//#####################################################################
// Constructor
//#####################################################################
LOG_ENTRY::
LOG_ENTRY(LOG_ENTRY* parent_input,const int depth_input,const int timer_id_input,
    const std::string& name_input,int& verbosity_level_input)
    :parent(parent_input),depth(depth_input),timer_id(timer_id_input),time(0),name(name_input),
    verbosity_level(verbosity_level_input)
{
    end_on_separate_line=false;log_file_end_on_separate_line=false;
    timer_start_time=TIMER::Singleton()->Peek_Time(timer_id);
}
//#####################################################################
// Destructor
//#####################################################################
LOG_ENTRY::
~LOG_ENTRY()
{
}
//#####################################################################
// Function Start
//#####################################################################
void LOG_ENTRY::
Start(LOG_CLASS& instance)
{
    if(depth<=verbosity_level){
        if(start_on_separate_line) ::std::putchar('\n');start_on_separate_line=needs_indent=true;
        ::std::printf("%*s%-*s",2*depth,"",50-2*depth,name.c_str());::std::fflush(stdout);}
    if(instance.log_file){
        if(log_file_start_on_separate_line) ::std::putc('\n',instance.log_file);log_file_start_on_separate_line=log_file_needs_indent=true;
        if(instance.xml) Start_XML(instance);
        else ::std::fprintf(instance.log_file,"%*s%-*s",2*depth,"",50-2*depth,name.c_str());
        ::std::fflush(instance.log_file);}
    timer_start_time=TIMER::Singleton()->Peek_Time(timer_id);
}
//#####################################################################
// Function Start_XML
//#####################################################################
void LOG_ENTRY::
Start_XML(LOG_CLASS& instance)
{
    ::std::fprintf(instance.log_file,"%*s<scope name=\"%s\">",2*depth,"",name.c_str());
}
//#####################################################################
// Function Stop
//#####################################################################
void LOG_ENTRY::
Stop(LOG_CLASS& instance)
{
    double time_since_start=(TIMER::Singleton()->Peek_Time(timer_id)-timer_start_time)/1000;
    if(depth<=verbosity_level){
        if(end_on_separate_line){
            if(start_on_separate_line) ::std::putchar('\n');
            ::std::printf("%*sEND %-*s",2*depth,"",50-2*depth-4,name.c_str());}
        end_on_separate_line=false;start_on_separate_line=needs_indent=true;
        ::std::printf("%8.4f s",time_since_start);::std::fflush(stdout);}
    if(instance.log_file){
        if(instance.xml){
            if(log_file_end_on_separate_line){
                if(log_file_start_on_separate_line) ::std::putc('\n',instance.log_file);
                ::std::fprintf(instance.log_file,"%*s",2*depth,"");}
            log_file_end_on_separate_line=false;log_file_start_on_separate_line=log_file_needs_indent=true;
            ::std::fprintf(instance.log_file,"<time value=\"%f\"/></scope>",time_since_start);::std::fflush(instance.log_file);}
        else{
            if(log_file_end_on_separate_line){
                if(log_file_start_on_separate_line) ::std::putc('\n',instance.log_file);
                ::std::fprintf(instance.log_file,"%*sEND %-*s",2*depth,"",50-2*depth-4,name.c_str());}
            log_file_end_on_separate_line=false;log_file_start_on_separate_line=log_file_needs_indent=true;
            ::std::fprintf(instance.log_file,"%8.4f s",time_since_start);::std::fflush(instance.log_file);}}
    time+=time_since_start;
}
//#####################################################################
// Function Get_Stop_Time
//#####################################################################
LOG_ENTRY* LOG_ENTRY::
Get_Stop_Time(LOG_CLASS& instance)
{
    Stop(instance);
    return parent;
}
//#####################################################################
// Function Get_New_Scope
//#####################################################################
LOG_ENTRY* LOG_ENTRY::
Get_New_Scope(LOG_CLASS& instance,const std::string& new_scope_identifier,const std::string& new_name)
{
    Stop(instance);
    return parent->Get_New_Scope(instance,new_scope_identifier,new_name);
}
//#####################################################################
// Function Get_New_Item
//#####################################################################
LOG_ENTRY* LOG_ENTRY::
Get_New_Item(LOG_CLASS& instance,const std::string& new_name)
{
    Stop(instance);
    return parent->Get_New_Item(instance,new_name);
}
//#####################################################################
// Function Get_Pop_Scope
//#####################################################################
LOG_ENTRY* LOG_ENTRY::
Get_Pop_Scope(LOG_CLASS& instance)
{
    Stop(instance);
    return parent->Get_Pop_Scope(instance);
}
//#####################################################################
// Function Dump_Log
//#####################################################################
void LOG_ENTRY::
Dump_Log(FILE* output)
{
    ::std::fprintf(output,"%*s%-*s%8.4f s\n",2*depth,"",50-2*depth,name.c_str(),time);
    ::std::fflush(output);
}
//#####################################################################
// Function Dump_Names
//#####################################################################
void LOG_ENTRY::
Dump_Names(FILE* output)
{
    ::std::fprintf(output,"%*s%-*s",2*depth,"",50-2*depth,name.c_str());
    ::std::fflush(output);
}
//#####################################################################
}
}
