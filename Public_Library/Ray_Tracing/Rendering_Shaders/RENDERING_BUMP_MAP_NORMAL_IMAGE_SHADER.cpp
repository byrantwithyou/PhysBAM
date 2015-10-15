//#####################################################################
// Copyright 2004, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Ray_Tracing/Rendering/PHOTON_MAP.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_SHADER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER<T>::
RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER(const MATERIAL_SHADER<T>* shader,RENDER_WORLD<T>& world,const T s0_input,
    const T s1_input,const T t0_input,const T t1_input)
    :RENDERING_BUMP_MAP_SHADER<T>(shader,world)
{
}
//#####################################################################
// Function Perturbed_Normal
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER<T>::
Perturbed_Normal(const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& object_space_position,
    const T s,const T t,const VECTOR<T,3>& object_space_normal,const VECTOR<T,3>& object_space_tangent,
    const VECTOR<T,3>& object_space_bitangent) const
{
    VECTOR<T,3> interpolated_normal=(interpolation.Clamped_To_Array_Cell(grid,normals,VECTOR<T,2>(s,t))).Normalized();
    MATRIX<T,3> tangent_space_to_object_space(object_space_tangent,object_space_bitangent,object_space_normal);
    return tangent_space_to_object_space*interpolated_normal;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER<T>::
Initialize(const std::string& filename,const T max_phi)
{
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > pixels;
    IMAGE<T>::Read(filename,pixels);grid.Initialize(pixels.Size(),RANGE<VECTOR<T,2> >::Unit_Box());
    normals.Resize(pixels.domain.min_corner.x,pixels.domain.max_corner.x,pixels.domain.min_corner.y,pixels.domain.max_corner.y);
    // convert each RGB value to its corresponding normal value
    for(int i=normals.domain.min_corner.x;i<normals.domain.max_corner.x;++i)
        for(int j=normals.domain.min_corner.y;j<normals.domain.max_corner.y;++j)
            normals(i,j)=pixels(i,j)*2-VECTOR<T,3>(1,1,1);
}
//#####################################################################
// Function Valid
//#####################################################################
template<class T> bool RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER<T>::
Valid() const
{
    return grid.counts.Min()>0 && grid.counts==normals.Size();
}
template class RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER<double>;
template class RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER<float>;
}
