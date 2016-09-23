//#####################################################################
// Copyright 2003-2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_TORRANCE_SPARROW_SHADER__
#define __RENDERING_TORRANCE_SPARROW_SHADER__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/constants.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Ray_Tracing/Rendering/FRESNEL.h>
#include <Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_TORRANCE_SPARROW_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    using MATERIAL_SHADER<T>::world;

    const MATERIAL_SHADER<T>& shader;

    T exponent,index_of_refraction,absorption_coefficient;

    RENDERING_TORRANCE_SPARROW_SHADER(const MATERIAL_SHADER<T>& shader_input,const T exponent_input,
        const T index_of_refraction_input,const T absorption_coefficient_input,RENDER_WORLD<T>& world_input);
    TV BRDF(const TV& reflected_direction,const TV& incident_direction,const TV& normal,const TV& reflectance) const;
    TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
        const TV& intersection_point,const TV& same_side_normal) const;
//#####################################################################
};
}
#endif
