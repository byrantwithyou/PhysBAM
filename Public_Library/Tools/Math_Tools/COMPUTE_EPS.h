//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPUTE_EPS
//#####################################################################
#ifndef __COMPUTE_EPS__
#define __COMPUTE_EPS__

namespace PhysBAM{
constexpr inline int compute_eps_helper(int r)
{return ((0x200880>>r)&1)-((0x088020>>r)&1);}

constexpr inline int compute_eps(int i,int j,int k)
{return compute_eps_helper(i+3*j+9*k);}
}
#endif
