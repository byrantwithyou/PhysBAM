//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLIP_LEVELSET__
#define __FLIP_LEVELSET__

#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>

namespace PhysBAM{

template<class TV>
class FLIP_LEVELSET_COMPLEMENT:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
public:
    const IMPLICIT_OBJECT<TV>& object;

    FLIP_LEVELSET_COMPLEMENT(const IMPLICIT_OBJECT<TV>& object)
        :object(object)
    {}

    ~FLIP_LEVELSET_COMPLEMENT()
    {delete &object;}

    T operator()(const TV& location) const
    {return -object(location);}

    T Extended_Phi(const TV& location) const
    {return -object(location);}

    TV Normal(const TV& location,const int aggregate=-1) const
    {return -object.Normal(location,aggregate);}
};

template<class TV>
class FLIP_LEVELSET_INTERSECTION:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
public:
    const IMPLICIT_OBJECT<TV>& object1;
    const IMPLICIT_OBJECT<TV>& object2;

    FLIP_LEVELSET_INTERSECTION(const IMPLICIT_OBJECT<TV>& object1,const IMPLICIT_OBJECT<TV>& object2)
        :object1(object1),object2(object2)
    {}

    ~FLIP_LEVELSET_INTERSECTION()
    {delete &object1;delete &object2;}

    T operator()(const TV& location) const
    {return max(object1(location),object2(location));}

    T Extended_Phi(const TV& location) const
    {return max(object1.Extended_Phi(location),object2.Extended_Phi(location));}

    TV Normal(const TV& location,const int aggregate=-1) const
    {if(object1(location)>object2(location)) return object1.Normal(location,aggregate);return object2.Normal(location,aggregate);}
};
}
#endif
