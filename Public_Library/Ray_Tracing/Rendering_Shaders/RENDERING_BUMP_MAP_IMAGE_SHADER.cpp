//#####################################################################
// Copyright 2004, Geoffrey Irving, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Tools/Images/IMAGE.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Ray_Tracing/Rendering/PHOTON_MAP.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_IMAGE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_SHADER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_BUMP_MAP_IMAGE_SHADER<T>::
RENDERING_BUMP_MAP_IMAGE_SHADER(const MATERIAL_SHADER<T>* shader,RENDER_WORLD<T>& world,const T s0_input,
    const T s1_input,const T t0_input,const T t1_input)
    :RENDERING_BUMP_MAP_SHADER<T>(shader,world),levelset(grid,phi),s0(s0_input),s_min(min(s0_input,s1_input)),
    s_max(max(s0_input,s1_input)),s_scaling_factor((T)1/(s1_input-s0_input)),t0(t0_input),
    t_min(min(t0_input,t1_input)),t_max(max(t0_input,t1_input)),t_scaling_factor((T)1/(t1_input-t0_input))
{
}
//#####################################################################
// Function Perturbed_Normal
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_BUMP_MAP_IMAGE_SHADER<T>::
Perturbed_Normal(const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& object_space_position,
    const T s,const T t,const VECTOR<T,3>& object_space_normal,const VECTOR<T,3>& object_space_tangent,
    const VECTOR<T,3>& object_space_bitangent) const
{
    if(s<s_min || s>s_max || t<t_min || t>t_max) return object_space_normal;
    VECTOR<T,2> location((s-s0)*s_scaling_factor,(t-t0)*t_scaling_factor);
    T fx=(levelset.Phi(VECTOR<T,2>(location.x+grid.dX.x,location.y))-levelset.Phi(VECTOR<T,2>(location.x-grid.dX.x,location.y)))/(2*grid.dX.x),
        fy=(levelset.Phi(VECTOR<T,2>(location.x,location.y+grid.dX.y))-levelset.Phi(VECTOR<T,2>(location.x,location.y-grid.dX.y)))/(2*grid.dX.y);
    VECTOR<T,3> N0(0,0,1),N1(-fx,-fy,1);
    return (object_space_normal+(N1-N0)).Normalized();
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> bool RENDERING_BUMP_MAP_IMAGE_SHADER<T>::
Initialize(const std::string& filename,const T max_phi)
{
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > pixels;
    IMAGE<T>::Read(filename,pixels);
    int i,j;//for(i=pixels.domain.min_corner.x;i<pixels.domain.max_corner.x;i++)for(j=pixels.domain.min_corner.y;j<pixels.domain.min_corner.y+pixels.domain.max_corner.y-j;j++) exchange(pixels(i,j),pixels(i,pixels.domain.max_corner.y+pixels.domain.min_corner.y-j));
    grid.Initialize(pixels.Size(),RANGE<VECTOR<T,2> >::Unit_Box()); 
    phi.Resize(pixels.domain);
    for(i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;++i)for(j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;++j)phi(i,j)=VECTOR<T,3>::Dot_Product(VECTOR<T,3>((T).299,(T).587,(T).114),pixels(i,j));
    phi*=max_phi/phi.Max();
    return true;
}
//#####################################################################
// Function Valid
//#####################################################################
template<class T> bool RENDERING_BUMP_MAP_IMAGE_SHADER<T>::
Valid() const
{
    return grid.counts.Min()>0 && grid.counts==phi.Size();
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void RENDERING_BUMP_MAP_IMAGE_SHADER<T>::
Print() const
{
    LOG::cout<<"grid "<<grid<<std::endl;
}
template class RENDERING_BUMP_MAP_IMAGE_SHADER<double>;
template class RENDERING_BUMP_MAP_IMAGE_SHADER<float>;
}
