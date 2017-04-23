//#####################################################################
// Copyright 2004-2008, Nipun Kwatra, Frank Losasso, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TIMER
//#####################################################################
#ifndef __TIMER__
#define __TIMER__

#include <Core/Arrays/ARRAY.h>
namespace PhysBAM{

double Get_Current_Time();
double Initialize_Timer();

class TIMER
{
private:
    struct DATA
    {
        double start,elapsed,accumulator;
    };

    double resolution;
    ARRAY<DATA> timers;
    ARRAY<int> free_timers;
    static TIMER* singleton_instance;
    double overhead;
public:

    TIMER();
    TIMER(const TIMER&) = delete;
    void operator=(const TIMER&) = delete;
    ~TIMER();

    //TODO: When to delete?
    static inline void Destroy_Singleton()
    {delete singleton_instance;singleton_instance=0;}

    static inline TIMER* Singleton()
    {if(!singleton_instance) singleton_instance=new TIMER;return singleton_instance;}

//#####################################################################
    double Get_Time();
    int Register_Timer();
    void Release_Timer(const int id);
    double Get_Total_Time_Since_Registration(const int id);
    double Peek_And_Reset_Time(const int id);
    void Reset_Time(const int id);
    double Peek_Time(const int id);
    void Start(const int id);
    void Stop(const int id);
    void Print_Stats(const int id,const char* str);
//#####################################################################
};
}
#endif
