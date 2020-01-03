//#####################################################################
// Copyright 2020, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Functions cdiv and fdiv
//#####################################################################
//
// guaranteed floor(a/b) and ceil(a/b), regardless of sign.
//               
//#####################################################################
#ifndef __int_div__
#define __int_div__

namespace PhysBAM{

// floor(a/b)
inline constexpr int fdiv(int a,int b)
{
    int d=a/b,r=a%b;
    return d-((r!=0)&((a^b)<0));
}

// ceil(a/b)
inline constexpr int cdiv(int a,int b)
{
    int d=a/b,r=a%b;
    return d+((r!=0)&((a^b)>=0));
}

// If b>0 then:

// b*c>=a implies c>=fdiv(a,b)
// c<fdiv(a,b) implies b*c<a
// b*c<=a equals c<=fdiv(a,b)
// b*c>a equals c>fdiv(a,b)

// b*c>=a equals c>=cdiv(a,b)
// b*c<a equals c<cdiv(a,b)
// c<=fdiv(a,b) implies b*c<=a
// b*c>a implies c>cdiv(a,b)

}
#endif

