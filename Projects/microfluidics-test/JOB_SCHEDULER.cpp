//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include "JOB_SCHEDULER.h"
#include <climits>
#include <fstream>
#include <mutex>
#include <thread>
#include <x86intrin.h>

namespace PhysBAM{
bool use_job_timing;

template<class T,class P> bool PRIORITY_QUEUE<T,P>::
Pop(T& data)
{
    std::unique_lock<std::mutex> lck(mtx);
    if(!a.m) return false;
    data=a(0).x;
    std::pop_heap(a.begin(),a.end(),cmp);
    a.Pop();
    return true;
}
template<class T,class P> int PRIORITY_QUEUE<T,P>::
Push(const T& data,P priority)
{
    std::unique_lock<std::mutex> lck(mtx);
    a.Append({data,priority});
    std::push_heap(a.begin(),a.end(),cmp);
    return a.m;
}

int THREAD_STATE::Get_Next_Job()
{
    int job=-1;
    if(priority_queue.Pop(job)) return job;

    while(1)
    {
        for(int i=0;i<core->threads.m;i++)
        {
            int t=(i+tid)%core->threads.m;
            if(core->threads(t)->priority_queue.Pop(job)) return job;
        }

        std::unique_lock<std::mutex> lck(core->mtx);
        if(++core->waiting==core->threads.m)
        {
            core->done.store(true);
            core->cv.notify_all();
            return -1;
        }
        core->cv.wait(lck);
        int w=--core->waiting;
        PHYSBAM_ASSERT(w>=0);
        if(core->done.load()) return -1;
    }
}

void THREAD_STATE::Release_Job(int job)
{
    int n=priority_queue.Push(job,core->jobs(job)->priority);
    if(n>1 && core->waiting.load())
        core->cv.notify_one();
}

void THREAD_STATE::Release_Dependencies(int job)
{
    for(auto a:core->jobs(job)->release_list)
        if(!--core->jobs(a)->missing_deps)
            Release_Job(a);
}

void Thread_Function(THREAD_STATE* ts)
{
    auto dur=[](auto a,auto b){return std::chrono::duration_cast<std::chrono::duration<double> >(b-a).count();};
    bool use_timing=use_job_timing;
    struct TIMING_DATA
    {
        double t0=0,t1=0;
        int job=0;
    };
    ARRAY<TIMING_DATA> timing_data;
    std::chrono::steady_clock::time_point t0,t1;
    if(use_timing) t0=std::chrono::steady_clock::now();
    
    while(1)
    {
        int job=ts->Get_Next_Job();
        if(job<0) break;

        std::chrono::steady_clock::time_point ta;
        if(use_timing) ta=std::chrono::steady_clock::now();

        ts->core->execute_job(ts->core->jobs(job)->job,ts->core->data);

        if(use_timing)
        {
            std::chrono::steady_clock::time_point tb=std::chrono::steady_clock::now();
            timing_data.Append({dur(t0,ta),dur(ta,tb),job});
            t0=tb;
        }

        ts->Release_Dependencies(job);
    }

    if(use_timing)
    {
        t1=std::chrono::steady_clock::now();
        std::ofstream fout(LOG::sprintf("timing-%i.txt",ts->tid));
        for(auto d:timing_data) LOG::fprintf(fout,"%i %.1f %.1f\n",d.job,d.t0*1e6,d.t1*1e6);
        LOG::fprintf(fout,"%.1f\n",dur(t0,t1)*1e6);
    }

}

void JOB_SCHEDULER_CORE::Execute_Jobs(int num_threads)
{
    if(num_threads<0) num_threads=std::thread::hardware_concurrency();
    LOG::printf("threads: %i\n",num_threads);
    LOG::printf("jobs: %i\n",jobs.m);

    threads.Resize(num_threads);
    for(int i=0;i<num_threads;i++)
    {
        threads(i)=new THREAD_STATE;
        threads(i)->core=this;
        threads(i)->tid=i;
    }

    for(int i=0,k=0;i<jobs.m;i++)
        if(!jobs(i)->missing_deps)
            threads(k++%num_threads)->Release_Job(i);

    for(auto& t:threads) t->th=new std::thread(Thread_Function,t);

    for(auto& t:threads)
    {
        t->th->join();
        delete t->th;
    }
}

int JOB_SCHEDULER_CORE::Compute_Priority_By_Paths_Rec(int j)
{
    if(jobs(j)->priority>=0)
        return jobs(j)->priority;
    int p=0;
    for(auto k:jobs(j)->release_list)
        p=std::max(p,Compute_Priority_By_Paths_Rec(k));
    p++;
    jobs(j)->priority=p;
    return p;
}

void JOB_SCHEDULER_CORE::Compute_Priority_By_Paths()
{
    for(auto& j:jobs) j->priority=-1;
     for(int i=jobs.m-1;i>=0;i--)
        Compute_Priority_By_Paths_Rec(i);
}

}
