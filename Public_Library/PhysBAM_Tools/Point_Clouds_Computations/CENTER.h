//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CENTER
//#####################################################################
#ifndef __CENTER__
#define __CENTER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
namespace PhysBAM{

namespace POINT_CLOUDS_COMPUTATIONS
{
    template<class TV,class T_ARRAY>
    static TV Center(const ARRAY_BASE<TV,T_ARRAY>& X)
    {typename TV::SCALAR total=X.Size();return total?X.Sum()/total:TV();}

    template<class TV,class T_ARRAY,class T_ARRAY2>
    static TV Weighted_Center(const ARRAY_BASE<TV,T_ARRAY>& X,const ARRAY_BASE<typename TV::SCALAR,T_ARRAY2>& weights)
    {typename TV::SCALAR total=weights.Sum();return total?X.Weighted_Sum(weights)/total:TV();}
}
}
#endif
