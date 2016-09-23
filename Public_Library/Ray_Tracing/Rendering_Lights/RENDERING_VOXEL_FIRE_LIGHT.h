//#####################################################################
// Copyright 2004-2007, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_VOXEL_FIRE_LIGHT
//#####################################################################
#ifndef __RENDERING_VOXEL_FIRE_LIGHT__
#define __RENDERING_VOXEL_FIRE_LIGHT__

#include <Core/Log/LOG.h>
#include <Core/Log/PROGRESS_INDICATOR.h>
#include <Core/Log/SCOPE.h>
#include <Core/Random_Numbers/PIECEWISE_CONSTANT_PDF.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_UNIFORM_VOXELS.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_VOXEL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_VOXEL_FIRE_LIGHT:public RENDERING_LIGHT<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
private:
    using RENDERING_LIGHT<T>::world;using RENDERING_LIGHT<T>::supports_global_photon_mapping;using RENDERING_LIGHT<T>::supports_caustic_photon_mapping;
    using RENDERING_LIGHT<T>::supports_volume_photon_mapping;using RENDERING_LIGHT<T>::global_photon_random;using RENDERING_LIGHT<T>::caustic_photon_random;
    using RENDERING_LIGHT<T>::volume_photon_random;using RENDERING_LIGHT<T>::sample_points_random;

    RENDERING_UNIFORM_VOXELS<T>& fire_voxels;
    RENDERING_VOXEL_SHADER<T>& fire_shader;
    T total_average_power;
    VECTOR<T,3> total_power;
    ARRAY<T,VECTOR<int,3> > pdf;
    ARRAY<PIECEWISE_CONSTANT_PDF<T> ,VECTOR<int,2> > z_cdf;
    ARRAY<PIECEWISE_CONSTANT_PDF<T> ,VECTOR<int,1> > y_cdf;
    PIECEWISE_CONSTANT_PDF<T> x_cdf;
public:

    RENDERING_VOXEL_FIRE_LIGHT(RENDERING_UNIFORM_VOXELS<T>& fire_voxels_input,
        RENDERING_VOXEL_SHADER<T>& fire_shader_input,RENDER_WORLD<T>& world_input,const bool supports_global_photons,
        const bool supports_caustic_photons,const bool supports_volume_photons);
    virtual ~RENDERING_VOXEL_FIRE_LIGHT();

    void Sample_Points(const VECTOR<T,3>& surface_position,const VECTOR<T,3>& surface_normal,
        ARRAY<RAY<VECTOR<T,3> > >& sample_array)const override;
    VECTOR<T,3> Sample_Point_In_Volume(const T xi_x,const T xi_y,const T xi_z,T& probability_of_location)const;
    VECTOR<T,3> Emitted_Light(const RENDERING_RAY<T>& ray) const override;
    int Emit_Photons(RENDERING_RAY<T>& parent_ray,PHOTON_MAP<T>& photon_map,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type) const  override;
//#####################################################################
};
}
#endif
