//#####################################################################
// Copyright 2017, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FINE_TIMER__
#define __FINE_TIMER__
#include <vector>
#include <x86intrin.h>

namespace PhysBAM
{
namespace FINE_TIMER
{
struct GLOBAL;
inline std::vector<GLOBAL*>& Get_Global()
{
    static std::vector<GLOBAL*> list;
    return list;
}

struct GLOBAL
{
    const char* name;
    long long total_time;
    long long child_time;
    long long num_calls;

    explicit GLOBAL(const char* name)
        :name(name),total_time(0),child_time(0),num_calls(0)
    {
        Get_Global().push_back(this);
    }
};

extern GLOBAL* current_scope;

struct LOCAL
{
    long long start_time;
    GLOBAL* global;
    GLOBAL* save_scope;

    explicit LOCAL(GLOBAL* global)
        :start_time(__rdtsc()),global(global),save_scope(current_scope)
    {
        current_scope=global;
    }

    ~LOCAL()
    {
        long long this_time=__rdtsc()-start_time;
        global->total_time+=this_time;
        global->num_calls++;
        if(save_scope) save_scope->child_time+=this_time;
        current_scope=save_scope;
    }
};

void Dump_Timing_Info();

#define TIMER_SCOPE_CONCAT_(a,b) a##b
#define TIMER_SCOPE_(n,c) static FINE_TIMER::GLOBAL TIMER_SCOPE_CONCAT_(global_,c)(n);FINE_TIMER::LOCAL TIMER_SCOPE_CONCAT_(local_,c)(&TIMER_SCOPE_CONCAT_(global_,c))
#define TIMER_SCOPE(n) TIMER_SCOPE_(n,__COUNTER__)
#define TIMER_SCOPE_FULL TIMER_SCOPE(__PRETTY_FUNCTION__)
#define TIMER_SCOPE_FUNC TIMER_SCOPE(__FUNCTION__)
}
}
#endif
