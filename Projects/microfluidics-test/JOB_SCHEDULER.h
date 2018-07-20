//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __JOB_SCHEDULER__
#define __JOB_SCHEDULER__
#include <Core/Arrays/ARRAY.h>
#include <mutex>
#include <condition_variable>

namespace PhysBAM{

struct JOB_SCHEDULER_CORE
{
    struct JOB_INFO
    {
        void* job;
        int priority;
        int missing_deps;
        ARRAY<int> release_list;
    };
    ARRAY<JOB_INFO> jobs;
    ARRAY<int> priority_queue;

    void (*execute_job)(void* job,void* data);
    void* data;

    std::mutex mtx;
    std::condition_variable cv;
    int jobs_left;
    
    int Add_Job(void* job,int priority)
    {return jobs.Append({job,priority,0});}

    // user depends on dep.
    void Register_Dependency(int dep,int user)
    {
        jobs(dep).release_list.Append(user);
        jobs(user).missing_deps++;
    }

    int Finish_Job(int job);
    void Execute_Jobs(int num_threads);
    void Compute_Priority_By_Paths();
    int Compute_Priority_By_Paths_Rec(int j);
};

template<class JOB,class DATA>
struct JOB_SCHEDULER:public JOB_SCHEDULER_CORE
{
    JOB_SCHEDULER(DATA* user_data)
    {
        execute_job=[](void* job,void* data){((JOB*)job)->Execute((DATA*)data);};
        data=user_data;
    }

    int Add_Job(JOB* job,int priority)
    {return JOB_SCHEDULER_CORE::Add_Job(job,priority);}
};

}
#endif
