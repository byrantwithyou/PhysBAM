//#####################################################################
// Copyright 2004, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER  
//#####################################################################
#ifndef __RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER__
#define __RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER__

#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER:public RENDERING_BUMP_MAP_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > normals;
    GRID<VECTOR<T,2> > grid;
    LINEAR_INTERPOLATION_UNIFORM<VECTOR<T,2>,VECTOR<T,3> > interpolation;

    RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER(const MATERIAL_SHADER<T>* shader,RENDER_WORLD<T>& world,
        const T s0_input=0,const T s1_input=1,const T t0_input=0,const T t1_input=1);

    virtual VECTOR<T,3> Perturbed_Normal(const RENDERING_OBJECT<T>& intersection_object,
        const VECTOR<T,3>& object_space_position,const T s,const T t,const VECTOR<T,3>& object_space_normal,
        const VECTOR<T,3>& object_space_tangent,const VECTOR<T,3>& object_space_bitangent) const;
    void Initialize(const std::string& filename,const T max_phi=.01);
    bool Valid() const;
};
}
#endif

