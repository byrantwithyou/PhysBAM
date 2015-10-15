//#####################################################################
// Copyright 2004, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_PHONG_SHADER__
#define __RENDERING_PHONG_SHADER__

#include <Tools/Math_Tools/max.h>
#include <Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_PHONG_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    using MATERIAL_SHADER<T>::world;

    TV ambient_color,diffuse_color,specular_color;
    T specular_exponent;
    const MATERIAL_SHADER<T>& shader;

    RENDERING_PHONG_SHADER(const TV& ambient_color_input,const TV& diffuse_color_input,
        const TV& specular_color_input,const T specular_exponent_input,
        const MATERIAL_SHADER<T>& shader_input,RENDER_WORLD<T>& world_input);
    RENDERING_PHONG_SHADER(const TV& ambient_color_input,const TV& diffuse_color_input,
        const T diffuse_coefficient,const TV& specular_color_input,const T specular_coefficient,
        const T specular_exponent_input,const MATERIAL_SHADER<T>& shader_input,RENDER_WORLD<T>& world_input);

    TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
        const TV& intersection_point,const TV& same_side_normal) const;
//#####################################################################
};
}
#endif
