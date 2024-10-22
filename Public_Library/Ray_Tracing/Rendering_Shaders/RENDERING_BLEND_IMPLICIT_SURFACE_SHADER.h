//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_BLEND_IMPLICIT_SURFACE_SHADER__
#define __RENDERING_BLEND_IMPLICIT_SURFACE_SHADER__

#include <Core/Math_Tools/INTERVAL.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BLEND_SHADER.h>
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;
template<class T>
class RENDERING_BLEND_IMPLICIT_SURFACE_SHADER:public RENDERING_BLEND_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    const IMPLICIT_OBJECT<TV>& implicit_surface;
    INTERVAL<T> value_range,weight_range;

    RENDERING_BLEND_IMPLICIT_SURFACE_SHADER(const IMPLICIT_OBJECT<TV>& implicit_surface_input,const INTERVAL<T>& value_range_input,const INTERVAL<T>& weight_range_input,
        const MATERIAL_SHADER<T>& shader1_input,const MATERIAL_SHADER<T>& shader2_input,RENDER_WORLD<T>& world_input)
        :RENDERING_BLEND_SHADER<T>(shader1_input,shader2_input,false,world_input),implicit_surface(implicit_surface_input),value_range(value_range_input),weight_range(weight_range_input)
    {}

    T Blending_Fraction(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {T phi=implicit_surface.Extended_Phi(intersection_point);
    return clamp((phi-value_range.min_corner)/value_range.Size(),weight_range.min_corner,weight_range.max_corner);}

//#####################################################################
};
}
#endif
