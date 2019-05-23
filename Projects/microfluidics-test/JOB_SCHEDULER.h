//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __JOB_SCHEDULER__
#define __JOB_SCHEDULER__
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace PhysBAM{

struct JOB_SCHEDULER_CORE;

struct alignas(64) JOB_INFO
{
    void* job;
    int priority;
    std::atomic<int> missing_deps;
    ARRAY<int> release_list;

    JOB_INFO(void* job,int priority)
        :job(job),priority(priority)
    {
        missing_deps.store(0);
    }
};

template<class T,class P>
struct PRIORITY_QUEUE
{
    ARRAY<PAIR<T,P> > a;
    std::mutex mtx;
    static bool cmp(const PAIR<T,P>& r,const PAIR<T,P>& s){return r.y<s.y;}
    
    bool Pop(T& data);
    int Push(const T& data,P priority);
};

struct THREAD_STATE
{
    JOB_SCHEDULER_CORE * core=0;

    PRIORITY_QUEUE<int,int> priority_queue;
    std::thread* th=0;
    int tid=-1;

    int Get_Next_Job();
    void Release_Dependencies(int job);
    void Release_Job(int job);
};

struct JOB_SCHEDULER_CORE
{
    ARRAY<THREAD_STATE*> threads;
    ARRAY<JOB_INFO*> jobs;

    std::mutex mtx;
    std::condition_variable cv;
    std::atomic<int> waiting;
    std::atomic<bool> done;

    void (*execute_job)(void* job,void* data);
    void* data;

    JOB_SCHEDULER_CORE()
    {
        waiting.store(0);
        done.store(false);
    }

    ~JOB_SCHEDULER_CORE()
    {
        jobs.Delete_Pointers_And_Clean_Memory();
        threads.Delete_Pointers_And_Clean_Memory();
    }
    
    int Add_Job(void* job,int priority)
    {return jobs.Append(new JOB_INFO(job,priority));}

    // user depends on dep.
    void Register_Dependency(int dep,int user)
    {
        jobs(dep)->release_list.Append(user);
        jobs(user)->missing_deps++;
    }

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
