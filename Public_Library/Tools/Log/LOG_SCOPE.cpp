//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Log/LOG_ENTRY.h>
#include <Tools/Log/LOG_SCOPE.h>
#include <Tools/Utilities/TIMER.h>
#include <string>
namespace PhysBAM{
namespace LOG{
LOG_SCOPE::
LOG_SCOPE(LOG_ENTRY* parent_input,int depth_input,int timer_id_input,const std::string& scope_identifier_input,
    const std::string& name_input,int& verbosity_level_input)
    :LOG_ENTRY(parent_input,depth_input,timer_id_input,name_input,verbosity_level_input),
    scope_identifier(scope_identifier_input)
{
}
LOG_SCOPE::
~LOG_SCOPE()
{
    children.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Get_Stop_Time
//#####################################################################
LOG_ENTRY* LOG_SCOPE::
Get_Stop_Time(LOG_CLASS& instance)
{
    return this;
}
//#####################################################################
// Function Get_New_Scope
//#####################################################################
LOG_ENTRY* LOG_SCOPE::
Get_New_Scope(LOG_CLASS& instance,const std::string& new_scope_identifier,const std::string& new_name)
{
    int entry;end_on_separate_line=true;log_file_end_on_separate_line=true;
    if(entries.Get(new_scope_identifier,entry)){
        children(entry)->name=new_name;
        return children(entry);}
    LOG_ENTRY* new_entry=new LOG_SCOPE(this,depth+1,timer_id,new_scope_identifier,new_name,verbosity_level);
    entries.Insert(new_scope_identifier,children.Append(new_entry));
    return new_entry;
}
//#####################################################################
// Function Get_New_Item
//#####################################################################
LOG_ENTRY* LOG_SCOPE::
Get_New_Item(LOG_CLASS& instance,const std::string& new_name)
{
    int entry;end_on_separate_line=true;log_file_end_on_separate_line=true;
    if(entries.Get(new_name,entry)) return children(entry);
    LOG_ENTRY* new_entry=new LOG_ENTRY(this,depth+1,timer_id,new_name,verbosity_level);
    entries.Insert(new_name,children.Append(new_entry));
    return new_entry;
}
//#####################################################################
// Function Get_Pop_Scope
//#####################################################################
LOG_ENTRY* LOG_SCOPE::
Get_Pop_Scope(LOG_CLASS& instance)
{
    Stop(instance);return parent;
}
//#####################################################################
// Function Start_XML
//#####################################################################
void LOG_SCOPE::
Start_XML(LOG_CLASS& instance)
{
    fprintf(instance.log_file,"%*s<scope id=\"%s\" name=\"%s\">",2*depth,"",scope_identifier.c_str(),name.c_str());
}
//#####################################################################
// Function Dump_Log
//#####################################################################
void LOG_SCOPE::
Dump_Log(FILE* output)
{
    fprintf(output,"%*s%-*s%8.4f\n",2*depth,"",50-2*depth,scope_identifier.c_str(),time);::std::fflush(output);
    for(int i=0;i<children.m;i++) children(i)->Dump_Log(output);
}
//#####################################################################
// Function Dump_Names
//#####################################################################
void LOG_SCOPE::
Dump_Names(FILE* output)
{
    LOG_ENTRY::Dump_Names(output);for(int i=0;i<children.m;i++)children(i)->Dump_Names(output);
}
}
}
