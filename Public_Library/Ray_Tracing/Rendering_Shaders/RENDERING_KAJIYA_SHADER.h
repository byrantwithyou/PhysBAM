//#####################################################################
// Copyright 2007, Joyce Pan.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_KAJIYA_SHADER__
#define __RENDERING_KAJIYA_SHADER__

#include <Core/Math_Tools/max.h>
#include <Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <cmath>
namespace PhysBAM{

template<class T>
class RENDERING_KAJIYA_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    using MATERIAL_SHADER<T>::world;

    TV diffuse_color,specular_color;
    T diffuse_coefficient,specular_coefficient,specular_exponent;
    
    RENDERING_KAJIYA_SHADER(const TV& diffuse_color_input,const T diffuse_coefficient_input,const TV& specular_color_input,
        const T specular_coefficient_input,const T specular_exponent_input,RENDER_WORLD<T>& world_input);

    TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
        const TV& intersection_point,const TV& same_side_normal) const;
//#####################################################################
};
}
#endif
