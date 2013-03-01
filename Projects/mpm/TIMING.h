//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TIMING
//#####################################################################
#ifndef __TIMING__
#define __TIMING__
#define TIMING_START ::PhysBAM::TIMING timing;timing.Start();
#define TIMING_END(x) timing.Print(x);
#include <PhysBAM_Tools/Log/LOG.h>
#include <iomanip>
#include <iostream>
#ifdef _WIN32
#include <windows.h>
#include <stdio.h>
namespace PhysBAM{
class TIMING
{
    LARGE_INTEGER m_time;
    LARGE_INTEGER m_freq;
public:
    TIMING()
    {
        QueryPerformanceFrequency(&m_freq);
    }

    void Start()
    {
        QueryPerformanceCounter(&m_time);
    }

    float Restart()
    {
        float val = Elapsed();
        QueryPerformanceCounter(&m_time);
        return val;
    }

    float Elapsed()
    {
        return now()-(*this);
    }

    void Print(const char* _name)
    {
        float time=Elapsed();
        LOG::cout<<_name<<" took "<<time<<" milliseconds"<<std::endl;
        fflush(stdout);
        Restart();
    }

    static inline TIMING Now() 
    {
        TIMING t;
        t.Start();
        return t;
    }

    float operator -(const TIMING& t1) 
    {
        return float((double)(m_time.QuadPart-t1.m_time.QuadPart)/(double)m_freq.QuadPart*1000.0);
    }
};
}
#else
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
namespace PhysBAM{
class TIMING
{
    timeval m_time;
public:
    void Start() 
    {
        gettimeofday(&m_time,NULL);
    }

    float Restart() 
    {
        float val = Elapsed();
        gettimeofday(&m_time,NULL);
        return val;
    }

    float Elapsed() 
    {
        return Now()-(*this);
    }

    void Print(const char* _name)
    {
        float time=Elapsed();
        LOG::cout<<std::setw(50)<<_name<<std::setw(6)<<" took "<<std::setw(20)<<time<<std::setw(18)<<" milliseconds"<<std::endl;
        fflush(stdout);
        Restart();
    }

    static inline TIMING Now() 
    {
        TIMING t;
        t.Start();
        return t;
    }

    float operator -(const TIMING& t1) 
    {
        return (float)1000.0f*(m_time.tv_sec-t1.m_time.tv_sec)+1.0e-3f*(m_time.tv_usec-t1.m_time.tv_usec);
    }
};
}
#endif
#endif 
