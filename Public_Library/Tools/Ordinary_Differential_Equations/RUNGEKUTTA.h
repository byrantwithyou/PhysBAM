//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Duc Nguyen, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RUNGEKUTTA
//#####################################################################
#ifndef __RUNGEKUTTA__
#define __RUNGEKUTTA__

namespace PhysBAM{

static const double rk_table[4][3][2]=
{
    {},
    {{0,1}},
    {{0,1},{.5,0}},
    {{0,1},{.75,-.5},{1./3,.5}}
};

template<class TV>
class RUNGEKUTTA
{
    typedef typename TV::SCALAR T;
public:
    int order; // 1, 2 or 3
    T dt; // size of time step
    T time; // current time
    int substep; // current substep
    TV& u; // variable being advanced in time
    TV u_copy; // a copy of the variable

    RUNGEKUTTA(TV& u,const int order,const T dt,const T time)
        :order(order),dt(dt),time(time),substep(0),u(u),u_copy(u)
    {}

    RUNGEKUTTA(const RUNGEKUTTA&) = delete;
    void operator=(const RUNGEKUTTA&) = delete;

    ~RUNGEKUTTA() {}

    void Next()
    {
        if(T a=(T)rk_table[order][substep][0]) u.Copy(a,u_copy,1-a,u);
        time+=(T)rk_table[order][substep++][1]*dt;
    }

    bool Valid()
    {return substep<order;}
//#####################################################################
};
}
#endif
