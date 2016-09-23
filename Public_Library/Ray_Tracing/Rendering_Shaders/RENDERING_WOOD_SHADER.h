//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// This is mainly from the book "Advanced Renderman"
//#####################################################################
#ifndef __RENDERING_WOOD_SHADER__
#define __RENDERING_WOOD_SHADER__

#include <Core/Random_Numbers/NOISE.h>
#include <Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_WOOD_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:

    VECTOR<T,3> color1,color2;

    T ring_frequency;
    T ring_noise,ring_noise_frequency;
    T trunk_wobble,trunk_wobble_frequency;
    T angular_wobble,angular_wobble_frequency;
    T grain_frequency;
    T grainy;
    T ringy;

    RENDERING_WOOD_SHADER(const VECTOR<T,3>& color1_input,const VECTOR<T,3>& color2_input,
        const T ring_frequency,const T ring_noise,const T ring_noise_frequency,const T trunk_wobble,
        const T trunk_wobble_frequency,const T angular_wobble,const T angular_wobble_frequency,
        const T grain_frequency,const T grainy,const T ringy,RENDER_WORLD<T>& world);

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,
        const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,
        const VECTOR<T,3>& same_side_normal) const;
//#####################################################################
};
}
#endif
