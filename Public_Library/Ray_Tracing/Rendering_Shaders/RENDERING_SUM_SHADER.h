//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_SUM_SHADER__
#define __RENDERING_SUM_SHADER__

#include <Core/Arrays/ARRAY.h>
#include <Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_SUM_SHADER:public MATERIAL_SHADER<T>
{
public:
    ARRAY<MATERIAL_SHADER<T>*> shaders;

    RENDERING_SUM_SHADER(RENDER_WORLD<T>& world_input);

    virtual VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,
        const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,
        const VECTOR<T,3>& same_side_normal) const;
    VECTOR<T,3> Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
        const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal,const RENDERING_LIGHT<T>& light,
        const RENDERING_RAY<T>& full_ray) const;
//#####################################################################
};
}
#endif
