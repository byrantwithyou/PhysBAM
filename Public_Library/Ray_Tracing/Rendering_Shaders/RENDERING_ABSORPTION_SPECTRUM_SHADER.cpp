//#####################################################################
// Copyright 2003-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_ABSORPTION_SPECTRUM_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BLEND_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
namespace PhysBAM{
template<class T> class RENDERING_BLEND_SHADER;

template<class T> RENDERING_ABSORPTION_SPECTRUM_SHADER<T>::
RENDERING_ABSORPTION_SPECTRUM_SHADER(const T absorption_coefficient_input,
    const VECTOR<T,3>& absorption_spectrum_input,RENDER_WORLD<T>& world_input,const T absorption_clamp_input)
    :VOLUMETRIC_SHADER<T>(world_input),absorption_coefficient1(absorption_coefficient_input),
    absorption_coefficient2(0),absorption_clamp(absorption_clamp_input),
    absorption_spectrum1(absorption_spectrum_input),blend_shader(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_ABSORPTION_SPECTRUM_SHADER<T>::
RENDERING_ABSORPTION_SPECTRUM_SHADER(const T absorption_coefficient1_input,const T absorption_coefficient2_input,
    const VECTOR<T,3>& absorption_spectrum1_input,const VECTOR<T,3>& absorption_spectrum2_input,
    const RENDERING_BLEND_SHADER<T>& blend_shader_input,RENDER_WORLD<T>& world_input,const T absorption_clamp_input)
    :VOLUMETRIC_SHADER<T>(world_input),absorption_coefficient1(absorption_coefficient1_input),
    absorption_coefficient2(absorption_coefficient2_input),absorption_clamp(absorption_clamp_input),
    absorption_spectrum1(absorption_spectrum1_input),absorption_spectrum2(absorption_spectrum2_input),
    blend_shader(&blend_shader_input)
{
}
//#####################################################################
// Function Compute_Absorption
//#####################################################################
template<class T> void RENDERING_ABSORPTION_SPECTRUM_SHADER<T>::
Compute_Absorption(T& absorption_coefficient,VECTOR<T,3>& absorption_spectrum,const RENDERING_RAY<T>& ray,
    const RENDERING_OBJECT<T>& object)
{
    if(!blend_shader){absorption_coefficient=absorption_coefficient1;absorption_spectrum=absorption_spectrum1;return;}
    VECTOR<T,3> sample_point=ray.ray.Point(ray.ray.t_max/2);
    // TODO: the following extreme hack assumes that the blender doesn't use most of it's parameters
    T weight=blend_shader->Blending_Fraction(ray,object,object,object,sample_point,VECTOR<T,3>());
    absorption_coefficient=(1-weight)*absorption_coefficient1+weight*absorption_coefficient2;
    absorption_spectrum=(1-weight)*absorption_spectrum1+weight*absorption_spectrum2;
}
//#####################################################################
// Function Attenuate_Color
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_ABSORPTION_SPECTRUM_SHADER<T>::
Attenuate_Color(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const VECTOR<T,3>& color)
{
    if(ray.ray.semi_infinite||!object.Inside(ray.ray.Point(ray.ray.t_max/2))) return color;
    T absorption_coefficient;
    VECTOR<T,3> absorption_spectrum;
    Compute_Absorption(absorption_coefficient,absorption_spectrum,ray,object);
    VECTOR<T,3> attenuation_coefficient=-ray.ray.t_max*absorption_coefficient*absorption_spectrum;
    for(int i=0;i<3;i++) if(absorption_clamp && abs(attenuation_coefficient(i))>absorption_clamp) attenuation_coefficient(i)=sign(attenuation_coefficient(i))*absorption_clamp;
    VECTOR<T,3> attenuation(exp(attenuation_coefficient.x),exp(attenuation_coefficient.y),exp(attenuation_coefficient.z));
    if(ray.debug_ray) ray.debug_ray->Add_Comment(LOG::sprintf("attenuation=%f %f %f\n",attenuation.x,attenuation.y,attenuation.z));
    return color*attenuation;
}
//#####################################################################
// Function Attenuate_Photon
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_ABSORPTION_SPECTRUM_SHADER<T>::
Attenuate_Photon(const RENDERING_RAY<T>& ray, const RENDERING_OBJECT<T>& object, const VECTOR<T,3>& photon_power, bool& should_throw)
{
    T absorption_coefficient;
    VECTOR<T,3> absorption_spectrum;
    Compute_Absorption(absorption_coefficient,absorption_spectrum,ray,object);
    VECTOR<T,3> attenuation_coefficient=-ray.ray.t_max*absorption_coefficient*absorption_spectrum;
    for(int i=0;i<3;i++) if(absorption_clamp && abs(attenuation_coefficient(i))>absorption_clamp) attenuation_coefficient(i)=sign(attenuation_coefficient(i))*absorption_clamp;
    VECTOR<T,3> attenuation(exp(attenuation_coefficient.x),exp(attenuation_coefficient.y),exp(attenuation_coefficient.z));
    if(ray.debug_ray) ray.debug_ray->Add_Comment(LOG::sprintf("photon attenuation=%f %f %f\n",attenuation.x,attenuation.y,attenuation.z));
    T average_probability=(T)((attenuation.x+attenuation.y+attenuation.z)/(T)3.0);
    T current_roll=world.random.Get_Uniform_Number((T)0.0,(T)1.0);
    if(current_roll<average_probability) should_throw=false;
    else should_throw=true;
    T scaling_factor=3/(attenuation.x+attenuation.y+attenuation.z);
    return photon_power*attenuation*scaling_factor;
}
//#####################################################################
// Function Attenuate_Light
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_ABSORPTION_SPECTRUM_SHADER<T>::
Attenuate_Light(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const RENDERING_LIGHT<T>& light,
    const VECTOR<T,3>& light_color)
{
    return Attenuate_Color(ray,object,light_color);
}
template class RENDERING_ABSORPTION_SPECTRUM_SHADER<double>;
template class RENDERING_ABSORPTION_SPECTRUM_SHADER<float>;
}
