//#####################################################################
// Copyright 2017, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/FINE_TIMER.h>
#include <Core/Log/LOG.h>
#include <algorithm>
#include <cstdio>
#include <cstring>

namespace PhysBAM
{
namespace FINE_TIMER
{
GLOBAL* current_scope;
void Dump_Timing_Info(bool sort_self)
{
    std::vector<GLOBAL*>& function_map=Get_Global();
    long long total_time=0;
    long long max_calls=0;
    for(const auto& it:function_map){
        total_time+=it->total_time-it->child_time;
        max_calls=std::max(max_calls,it->num_calls);}

    std::sort(function_map.begin(),function_map.end(),
        [sort_self](const GLOBAL*a,const GLOBAL*b)
        {
            if(sort_self){
                long long as=a->total_time-a->child_time;
                long long bs=b->total_time-b->child_time;
                if(as!=bs) return as>bs;}
            if(a->total_time!=b->total_time) return a->total_time>b->total_time;
            if(a->num_calls!=b->num_calls) return a->num_calls>b->num_calls;
            return strcmp(a->name,b->name)<0;
        });

    int calls_digits=1;
    while(max_calls/=10) calls_digits++;
    calls_digits=std::max(calls_digits,5);
    LOG::printf("\ntotal   self  %*s name\n",calls_digits,"calls");
    for(const auto& it:function_map)
    {
        LOG::printf("%.4f %.4f %*lli %s\n",
            (double)it->total_time/total_time,
            (double)(it->total_time-it->child_time)/total_time,
            calls_digits,it->num_calls,it->name);
    }
}
void Dump_Timing_Info()
{
    Dump_Timing_Info(false);
    Dump_Timing_Info(true);
}
}
}
