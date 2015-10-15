//#####################################################################
// Copyright 2003-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Ray_Tracing/Rendering/PHOTON_MAP.h>
#include <Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_COLOR_BLEND_SHADER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_COLOR_BLEND_SHADER<T>::
RENDERING_COLOR_BLEND_SHADER(RENDER_WORLD<T>& world_input) 
    :MATERIAL_SHADER<T>(world_input),total_weight(0)
{
}
//#####################################################################
// Function Add_Shader
//#####################################################################
template<class T> void RENDERING_COLOR_BLEND_SHADER<T>::
Add_Shader(const MATERIAL_SHADER<T>* shader,T weight)
{
    shaders.Append(shader);
    weights.Append(weight);
    total_weight+=weight;
}
//#####################################################################
// Function Evaluate_Diffuse_BRDF
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_COLOR_BLEND_SHADER<T>::
Evaluate_Diffuse_BRDF(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
    const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
    const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
{
    PHYSBAM_FATAL_ERROR();
#if 0
    VECTOR<T,3> color_accumulator;
    for(int i=0;i<shaders.m;i++)
        color_accumulator+=weights(i)*shaders(i)->Evaluate_Diffuse_BRDF(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return color_accumulator;
#endif
    return VECTOR<T,3>();
}
//#####################################################################
// Function Shade_Surface_Using_Direct_Illumination
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_COLOR_BLEND_SHADER<T>::
Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
    const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
    const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
{
    VECTOR<T,3> color_accumulator;
    for(int i=0;i<shaders.m;i++)color_accumulator+=weights(i)*shaders(i)->Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal); 
    return color_accumulator;
}
//#####################################################################
// Function Shade_Light_Ray
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_COLOR_BLEND_SHADER<T>::
Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
    const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
    const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,
    const RENDERING_RAY<T>& full_ray) const
{
    VECTOR<T,3> color_accumulator;
    for(int i=0;i<shaders.m;i++)color_accumulator+=weights(i)*shaders(i)->Shade_Light_Ray(ray,exiting_object,entering_object,intersection_object,intersection_point,
        same_side_normal,light,full_ray);
    return color_accumulator;
}
//#####################################################################
// Function Shade_Surface_Using_Indirect_Illumination
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_COLOR_BLEND_SHADER<T>::
Shade_Surface_Using_Indirect_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
    const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
    const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
{
    VECTOR<T,3> color_accumulator;
    for(int i=0;i<shaders.m;i++)color_accumulator+=weights(i)*shaders(i)->Shade_Surface_Using_Indirect_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal); 
    return color_accumulator;
}
//#####################################################################
// Function Shade_Surface_Using_Approximate_Full_Illumination
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_COLOR_BLEND_SHADER<T>::
Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,
    const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
    const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,
    const VECTOR<T,3>& same_side_normal) const
{
    VECTOR<T,3> color_accumulator;
    for(int i=0;i<shaders.m;i++)color_accumulator+=weights(i)*shaders(i)->Shade_Surface_Using_Approximate_Full_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    return color_accumulator;
}
//#####################################################################
// Function Receive_Photon
//#####################################################################
template<class T> void RENDERING_COLOR_BLEND_SHADER<T>::
Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
    const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
    const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const VECTOR<T,3>& power,
    const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const
{
    T r=world.random.Get_Uniform_Number((T)0,total_weight);
    for(int i=0;i<shaders.m-1;i++){
        if(r<=weights(i)){
            shaders(i)->Receive_Photon(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,power,type,diffuse_bounces,specular_bounces);
            return;}
        else r-=weights(i);}
    shaders.Last()->Receive_Photon(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal,power,type,diffuse_bounces,specular_bounces);
}
template class RENDERING_COLOR_BLEND_SHADER<double>;
template class RENDERING_COLOR_BLEND_SHADER<float>;
}
