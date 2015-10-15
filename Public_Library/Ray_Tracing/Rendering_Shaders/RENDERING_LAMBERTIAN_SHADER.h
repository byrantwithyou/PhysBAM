//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_LAMBERTIAN_SHADER__
#define __RENDERING_LAMBERTIAN_SHADER__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Math_Tools/constants.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_LAMBERTIAN_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    using MATERIAL_SHADER<T>::world;

    const MATERIAL_SHADER<T>& shader;
    T ambient_coefficient,diffuse_coefficient;
    TV ambient_color,diffuse_color;
    TV ambient_factor,diffuse_factor; // precomputed multiplications
    bool visualize_photon_map_directly;

    RENDERING_LAMBERTIAN_SHADER(const MATERIAL_SHADER<T>& shader_input,RENDER_WORLD<T>& world_input,
        bool visualize_photon_map_directly_input=false);
    RENDERING_LAMBERTIAN_SHADER(const MATERIAL_SHADER<T>& shader_input,const T ambient_coefficient_input,
        const TV& ambient_color_input,const T diffuse_coefficient_input,const TV& diffuse_color_input,
        RENDER_WORLD<T>& world_input,bool visualize_photon_map_directly_input=false);
    TV Evaluate_Diffuse_BRDF(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
        const TV& intersection_point,const TV& same_side_normal) const;
    TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
        const TV& intersection_point,const TV& same_side_normal) const;
    TV Shade_Surface_Using_Indirect_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
        const TV& intersection_point,const TV& same_side_normal) const;
    TV Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,
        const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const;
    void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,
        const RENDERING_OBJECT<T>& entering_object,const RENDERING_OBJECT<T>& intersection_object,
        const TV& intersection_point,const TV& same_side_normal,const TV& power,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const;
//#####################################################################
};
}
#endif
