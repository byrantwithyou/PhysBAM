//#####################################################################
// Copyright 2004-2008, Geoffrey Irving, Frank Losasso, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOG
//##################################################################### 
#ifndef __LOG__
#define __LOG__

#include <Tools/Utilities/NONCOPYABLE.h>
#include <cassert>
#include <ostream>
#include <sstream>
namespace PhysBAM{

class TIMER;

namespace LOG{

class LOG_ENTRY;
class LOG_SCOPE;

class LOG_CLASS
{
    friend class LOG_ENTRY;
    friend class LOG_SCOPE;
    friend class LOG_COUT_BUFFER;
    friend class LOG_CERR_BUFFER;
    friend void Reset();
    friend void Dump_Log();

    TIMER* timer_singleton;
    int timer_id;
    bool suppress_cout;
    bool suppress_cerr;
public:
    bool suppress_timing;
    FILE* log_file;
    int verbosity_level;
    bool log_file_temporary;
    bool xml;

    LOG_ENTRY* root;
    LOG_ENTRY* current_entry;

    LOG_CLASS(const bool suppress_cout,const bool suppress_cerr,const bool suppress_timing,const int verbosity_level,const bool cache_initial_output);
    ~LOG_CLASS();

    static void Push_Scope(const std::string& scope_identifier,const std::string& scope_name);
    static void Pop_Scope();
public:
    static void Time_Helper(const std::string& label);
    void Copy_Log_To_File(const std::string& filename,const bool append);

//##################################################################### 
};

// These next few lines are important to ensure no static data from LOG.cpp is accessed for DLLs
LOG_CLASS* Instance();
std::ostream& cout_Helper();
std::ostream& cerr_Helper();
extern std::ostream& cout;
extern std::ostream& cerr;

void Initialize_Logging(const bool suppress_cout_input=false,const bool suppress_timing_input=false,const int verbosity_level_input=1<<30,const bool cache_initial_output=false);
void Finish_Logging();
void Stop_Time();

//#####################################################################
// Function Stat
//#####################################################################
void Stat_Helper(const std::string& label,const std::stringstream& s);
template<class T_VALUE> void
Stat(const std::string& label,const T_VALUE& value)
{std::stringstream s;s<<value;Stat_Helper(label,s);}
void Reset(); 
void Dump_Log();

inline void Time(const std::string& format)
{if(Instance()->suppress_timing) return;
LOG_CLASS::Time_Helper(format);}

//##################################################################### 
}
}
#include <Tools/Log/LOG_PRINTF.h>
#endif
