//#####################################################################
// Copyright 2005, Jiayi Chong, Jeong-Mo Hong, Frank Losasso, Andy Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT
//#####################################################################
#ifndef __RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT__
#define __RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT__

#include <Core/Math_Tools/RANGE.h>
#include <Incompressible/Level_Sets/LEVELSET_MULTIPLE.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T,class T_LEVELSET_MULTIPLE>
class RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
    using RENDERING_OBJECT<T>::Inside; // silence -Woverloaded-virtual
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::inverse_transform;using RENDERING_OBJECT<T>::priority;

    T_LEVELSET_MULTIPLE& levelset_multiple;
    int region;
    
    RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT(T_LEVELSET_MULTIPLE& levelset_multiple_input, int region_input)
        :levelset_multiple(levelset_multiple_input),region(region_input)
    {}

    virtual ~RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT()
    {}

    TV Normal(const TV& location,const int aggregate=0) const override
    {assert((aggregate >= 1 && aggregate <= 6) || aggregate == -1);
    if(aggregate != -1) return RANGE<TV>(levelset_multiple.grid.domain).Normal(aggregate);else return levelset_multiple.levelsets(region)->Normal(location);}

    bool Inside(const TV& location) const override
    {if(!levelset_multiple.grid.domain.Inside(location,small_number)) return false;
    return region==levelset_multiple.Inside_Region(location);}

//#####################################################################
};   
}
#endif

