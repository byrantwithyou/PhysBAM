//#####################################################################
// Copyright 2003-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_ABSORPTION_SPECTRUM_SHADER__
#define __RENDERING_ABSORPTION_SPECTRUM_SHADER__

#include <Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
namespace PhysBAM{
template<class T> class RENDERING_BLEND_SHADER;

template<class T>
class RENDERING_ABSORPTION_SPECTRUM_SHADER:public VOLUMETRIC_SHADER<T>
{
public:
    using VOLUMETRIC_SHADER<T>::world;

    T absorption_coefficient1,absorption_coefficient2;
    T absorption_clamp;
    const VECTOR<T,3> absorption_spectrum1,absorption_spectrum2;
    const RENDERING_BLEND_SHADER<T>* blend_shader;
    
    RENDERING_ABSORPTION_SPECTRUM_SHADER(const T absorption_coefficient_input,
        const VECTOR<T,3>& absorption_spectrum_input,RENDER_WORLD<T>& world_input,const T absorption_clamp_input=0);

    RENDERING_ABSORPTION_SPECTRUM_SHADER(const T absorption_coefficient1_input,const T absorption_coefficient2_input,
        const VECTOR<T,3>& absorption_spectrum1_input,const VECTOR<T,3>& absorption_spectrum2_input,
        const RENDERING_BLEND_SHADER<T>& blend_shader_input,RENDER_WORLD<T>& world_input,const T absorption_clamp_input=0);

    void Compute_Absorption(T& absorption_coefficient,VECTOR<T,3>& absorption_spectrum,const RENDERING_RAY<T>& ray,
        const RENDERING_OBJECT<T>& object);
    VECTOR<T,3> Attenuate_Color(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,
        const VECTOR<T,3>& color) override;
    VECTOR<T,3> Attenuate_Photon(const RENDERING_RAY<T>& ray, const RENDERING_OBJECT<T>& object,
        const VECTOR<T,3>& photon_power, bool& should_throw) override;
    VECTOR<T,3> Attenuate_Light(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,
        const RENDERING_LIGHT<T>& light,const VECTOR<T,3>& light_color) override;
//#####################################################################
};
}
#endif
