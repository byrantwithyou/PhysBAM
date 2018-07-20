//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <climits>
#include "JOB_SCHEDULER.h"
#include <mutex>
#include <thread>

namespace PhysBAM{

void Thread_Function(JOB_SCHEDULER_CORE* s)
{
    int job=-1;
    while(1)
    {
        job=s->Finish_Job(job);
        if(job<0) break;
        s->execute_job(s->jobs(job).job,s->data);
    }
}

int JOB_SCHEDULER_CORE::Finish_Job(int job)
{
    std::unique_lock<std::mutex> lck(mtx);
    if(job>=0){
        if(!--jobs_left){
            cv.notify_all();
            return -1;}
        int best_priority=INT_MIN;
        int best_job=-1;
        if(priority_queue.m){
            best_job=-1;
            best_priority=jobs(priority_queue(0)).priority;}

        for(auto a:jobs(job).release_list){
            if(--jobs(a).missing_deps>0) continue;
            if(jobs(a).priority>=best_priority){
                best_priority=jobs(a).priority;
                std::swap(best_job,a);
                if(a<0) continue;}
            priority_queue.Append(a);
            std::push_heap(priority_queue.begin(),priority_queue.end(),
                [this](int r,int s){return jobs(r).priority<jobs(s).priority;});
            cv.notify_one();}

        if(best_job>=0) return best_job;}

    while(!priority_queue.m && jobs_left) cv.wait(lck);
    if(!jobs_left) return -1;
    int j=priority_queue(0);
    std::pop_heap(priority_queue.begin(),priority_queue.end(),[this](int r,int s){return jobs(r).priority<jobs(s).priority;});
    priority_queue.Pop();
    return j;
}

void JOB_SCHEDULER_CORE::Execute_Jobs(int num_threads)
{
    if(num_threads<0) num_threads=std::thread::hardware_concurrency();
    printf("threads: %i\n",num_threads);
    jobs_left=jobs.m;
    for(int i=0;i<jobs.m;i++)
        if(!jobs(i).missing_deps)
            priority_queue.Append(i);
    
    std::make_heap(priority_queue.begin(),priority_queue.end(),
        [this](int r,int s){return jobs(r).priority<jobs(s).priority;});

    ARRAY<std::thread*> threads;
    for(int i=0;i<num_threads;i++) threads.Append(new std::thread(Thread_Function,this));
    
    for(int i=0;i<threads.m;i++){
        threads(i)->join();
        delete threads(i);}
}

int JOB_SCHEDULER_CORE::Compute_Priority_By_Paths_Rec(int j)
{
    if(jobs(j).priority>=0)
        return jobs(j).priority;
    int p=0;
    for(auto k:jobs(j).release_list)
        p=std::max(p,Compute_Priority_By_Paths_Rec(k));
    jobs(j).priority=p;
    return p;
}

void JOB_SCHEDULER_CORE::Compute_Priority_By_Paths()
{
    for(auto& j:jobs) j.priority=-1;
     for(int i=jobs.m-1;i>=0;i--)
        Compute_Priority_By_Paths_Rec(i);
}

// //#####################################################################
// // Function Queue_Users
// //#####################################################################
// template<class T> void CACHED_ELIMINATION_MATRIX<T>::
// Queue_Users(int j)
// {
// #pragma omp critical
// }
// //#####################################################################
// // Function Get_Job
// //#####################################################################
// template<class T> int CACHED_ELIMINATION_MATRIX<T>::
// Get_Job()
// {
// #pragma omp critical
//     int j=priority_queue(0);
//     std::pop_heap(priority_queue.begin(),priority_queue.end(),[this](int r,int s){return jobs(r).priority<jobs(s).priority;});
//     priority_queue.Pop();
//     return j;
// }



}
