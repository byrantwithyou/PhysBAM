//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOG_SCOPE
//##################################################################### 
#ifndef __LOG_SCOPE__
#define __LOG_SCOPE__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Log/LOG_ENTRY.h>
#include <string>
namespace PhysBAM{
namespace LOG{

class LOG_SCOPE:public LOG_ENTRY
{
public:
    HASHTABLE<std::string,int> entries;
    ARRAY<LOG_ENTRY*> children;
    std::string scope_identifier;

    LOG_SCOPE(LOG_ENTRY* parent_input,int depth_input,int timer_id_input,const std::string& scope_identifier_input,
        const std::string& name_input,int& verbosity_level_input);
    virtual ~LOG_SCOPE();

    virtual LOG_ENTRY* Get_Stop_Time(LOG_CLASS& instance);
    LOG_ENTRY* Get_New_Scope(LOG_CLASS& instance,const std::string& new_scope_identifier,const std::string& new_name);
    LOG_ENTRY* Get_New_Item(LOG_CLASS& instance,const std::string& new_name);
    LOG_ENTRY* Get_Pop_Scope(LOG_CLASS& instance);
    void Start_XML(LOG_CLASS& instance) override;
    void Dump_Log(FILE* output) override;
    void Dump_Names(FILE* output) override;

//##################################################################### 
};
}
}
#endif
