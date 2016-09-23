//#####################################################################
// Copyright 2004, Geoffrey Irving, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_BUMP_MAP_IMAGE_SHADER  
//#####################################################################
#ifndef __RENDERING_BUMP_MAP_IMAGE_SHADER__
#define __RENDERING_BUMP_MAP_IMAGE_SHADER__

#include <Core/Log/LOG.h>
#include <Tools/Images/IMAGE.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_BUMP_MAP_IMAGE_SHADER:public RENDERING_BUMP_MAP_SHADER<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<T,2> TV2;
public:
    GRID<TV2> grid;
    ARRAY<T,VECTOR<int,2> > phi;
    LEVELSET<TV2> levelset;
    T s0,s_min,s_max,s_scaling_factor,t0,t_min,t_max,t_scaling_factor;

    RENDERING_BUMP_MAP_IMAGE_SHADER(const MATERIAL_SHADER<T>* shader,RENDER_WORLD<T>& world,const T s0_input=0,
        const T s1_input=1,const T t0_input=0,const T t1_input=1);
    VECTOR<T,3> Perturbed_Normal(const RENDERING_OBJECT<T>& intersection_object,
        const VECTOR<T,3>& object_space_position,const T s,const T t,const VECTOR<T,3>& object_space_normal,
        const VECTOR<T,3>& object_space_tangent,const VECTOR<T,3>& object_space_bitangent) const override;
    bool Initialize(const std::string& filename,const T max_phi=.01);
    bool Valid() const;
    void Print() const;
//#####################################################################
};
}
#endif

