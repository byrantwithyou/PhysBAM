//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Ray_Tracing/Rendering/PHOTON_MAP.h>
#include <Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_SUM_SHADER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_SUM_SHADER<T>::
RENDERING_SUM_SHADER(RENDER_WORLD<T>& world_input)
    :MATERIAL_SHADER<T>(world_input)
{
}
//#####################################################################
// Function Shade_Surface_Using_Direct_Illumination
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_SUM_SHADER<T>::
Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
    const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
    const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
{
    VECTOR<T,3> total_radiance;
    for(int i=0;i<shaders.m;i++) 
        total_radiance+=shaders(i)->Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return total_radiance;
}
//#####################################################################
// Function Shade_Light_Ray
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_SUM_SHADER<T>::
Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
    const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
    const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,
    const RENDERING_RAY<T>& full_ray) const
{
    VECTOR<T,3> total_radiance;
    for(int i=0;i<shaders.m;i++)
        total_radiance+=shaders(i)->Shade_Light_Ray(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,light,full_ray);
    return total_radiance;
}
template class RENDERING_SUM_SHADER<double>;
template class RENDERING_SUM_SHADER<float>;
}
