//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __MPLE_DOUBLE_WELL__
#define __MPLE_DOUBLE_WELL__

namespace PhysBAM{

template<class T>
class MPLE_DOUBLE_WELL
{
public:

    static T Gradient(const T x)
    {return  2*x*(x-1)*(2*x-1);}
};
}
#endif
