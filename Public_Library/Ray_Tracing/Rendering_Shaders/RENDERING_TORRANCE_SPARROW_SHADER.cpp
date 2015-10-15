//#####################################################################
// Copyright 2003-2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Math_Tools/constants.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Ray_Tracing/Rendering/FRESNEL.h>
#include <Ray_Tracing/Rendering/PHOTON_MAP.h>
#include <Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_TORRANCE_SPARROW_SHADER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_TORRANCE_SPARROW_SHADER<T>::
RENDERING_TORRANCE_SPARROW_SHADER(const MATERIAL_SHADER<T>& shader_input,const T exponent_input,
    const T index_of_refraction_input,const T absorption_coefficient_input,RENDER_WORLD<T>& world_input)
    :MATERIAL_SHADER<T>(world_input),shader(shader_input),exponent(exponent_input),
    index_of_refraction(index_of_refraction_input),absorption_coefficient(absorption_coefficient_input)
{}
//#####################################################################
// Function BRDF
//#####################################################################
template<class T> auto RENDERING_TORRANCE_SPARROW_SHADER<T>::
BRDF(const TV& reflected_direction,const TV& incident_direction,const TV& normal,const TV& reflectance) const -> TV
{
    T cos_reflected_normal=abs(TV::Dot_Product(normal,reflected_direction)),cos_incident_normal=abs(TV::Dot_Product(normal,incident_direction));
    TV half=(T).5*(reflected_direction+incident_direction);T cos_half_reflected=abs(TV::Dot_Product(half,reflected_direction));
    T cos_half_normal=abs(TV::Dot_Product(half,normal));
    T geometry_attenuation=min((T)1,min((T)2*cos_half_normal*cos_incident_normal/cos_half_reflected,(T)2*cos_half_normal*cos_reflected_normal/cos_half_reflected));
    T distribution=(exponent+2)*(1/(2*(T)pi))*(T)pow(max((T)0,cos_half_normal),exponent);
    T fresnel=FRESNEL<T>::Conductor_Reflection(index_of_refraction,cos_half_normal,absorption_coefficient);
    return fresnel*geometry_attenuation*distribution*reflectance/(4*cos_reflected_normal*cos_incident_normal);}
//#####################################################################
// Function Shade_Surface_Using_Direct_Illumination
//#####################################################################
template<class T> auto RENDERING_TORRANCE_SPARROW_SHADER<T>::
Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
    const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
    const TV& intersection_point,const TV& same_side_normal) const -> TV
{
    const ARRAY<RENDERING_LIGHT<T> *>& lights=world.Lights();
    TV same_side_position=intersection_point+same_side_normal*intersection_object.small_number*2;
    TV accumulated_color(0,0,0);
    TV reflectance=shader.Shade_Surface_Using_Direct_Illumination(ray,exiting_object,entering_object,intersection_object,intersection_point,same_side_normal);
    for(int light_index=0;light_index<lights.m;light_index++){
        TV accumulated_samples(0,0,0);
        ARRAY<RAY<TV> > sample_array;
        lights(light_index)->Sample_Points(same_side_position,same_side_normal,sample_array);
        for(int sample=0;sample<sample_array.m;sample++){
            RENDERING_RAY<T> ray_to_light(sample_array(sample),1,&exiting_object);
            TV light_color=world.Incident_Light(ray_to_light,*lights(light_index),ray_to_light,ray);
            T L_N=TV::Dot_Product(ray_to_light.ray.direction,same_side_normal);
            if(L_N<0) continue;
            accumulated_samples+=L_N*light_color*BRDF(-ray.ray.direction,ray_to_light.ray.direction,same_side_normal,reflectance);}
        accumulated_color+=accumulated_samples/T(sample_array.m);}
    return accumulated_color;
}
template class RENDERING_TORRANCE_SPARROW_SHADER<double>;
template class RENDERING_TORRANCE_SPARROW_SHADER<float>;
}
