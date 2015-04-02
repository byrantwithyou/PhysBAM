//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOG_ENTRY
//##################################################################### 
#ifndef __LOG_ENTRY__
#define __LOG_ENTRY__

#include <Tools/Log/LOG.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <cstdio>
#include <string>
namespace PhysBAM{
namespace LOG{

class LOG_ENTRY:public NONCOPYABLE
{
public:
    LOG_ENTRY* parent;
    int depth;
    int timer_id;
    double time;
    double timer_start_time;
    std::string name;
    bool end_on_separate_line,log_file_end_on_separate_line;
    static bool start_on_separate_line,log_file_start_on_separate_line;
    static bool needs_indent,log_file_needs_indent;
    int& verbosity_level;

    LOG_ENTRY(LOG_ENTRY* parent_input,const int depth_input,const int timer_id_input,
        const std::string& name_input,int& verbosity_level_input);
    virtual ~LOG_ENTRY();

    void Start(LOG_CLASS& instance);
    virtual void Start_XML(LOG_CLASS& instance);
    void Stop(LOG_CLASS& instance);
    virtual LOG_ENTRY* Get_Stop_Time(LOG_CLASS& instance);
    virtual LOG_ENTRY* Get_New_Scope(LOG_CLASS& instance,const std::string& new_scope_identifier,const std::string& new_name);
    virtual LOG_ENTRY* Get_New_Item(LOG_CLASS& instance,const std::string& new_name);
    virtual LOG_ENTRY* Get_Pop_Scope(LOG_CLASS& instance);
    virtual void Dump_Log(FILE* output);
    virtual void Dump_Names(FILE* output);

//##################################################################### 
};
}
}
#endif
