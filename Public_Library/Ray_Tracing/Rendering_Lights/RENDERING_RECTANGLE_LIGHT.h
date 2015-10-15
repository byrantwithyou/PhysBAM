//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_RECTANGLE_LIGHT
//#####################################################################
#ifndef __RENDERING_RECTANGLE_LIGHT__
#define __RENDERING_RECTANGLE_LIGHT__

#include <Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class T>
class RENDERING_RECTANGLE_LIGHT:public RENDERING_LIGHT<T>,public RENDERING_TRIANGULATED_SURFACE<T>
{
private:
    using RENDERING_LIGHT<T>::world;using RENDERING_LIGHT<T>::position;using RENDERING_LIGHT<T>::supports_global_photon_mapping;using RENDERING_LIGHT<T>::supports_caustic_photon_mapping;
    using RENDERING_LIGHT<T>::supports_volume_photon_mapping;using RENDERING_LIGHT<T>::color;using RENDERING_LIGHT<T>::brightness;using RENDERING_TRIANGULATED_SURFACE<T>::small_number;
    using RENDERING_TRIANGULATED_SURFACE<T>::triangulated_surface;using RENDERING_LIGHT<T>::global_photon_random;using RENDERING_LIGHT<T>::caustic_photon_random;
    using RENDERING_LIGHT<T>::volume_photon_random;using RENDERING_LIGHT<T>::sample_points_random;

    VECTOR<T,3> u_direction,v_direction;
    VECTOR<T,3> normal;
    const int u_samples,v_samples;
    T one_over_u_samples,one_over_v_samples;
    bool use_stratified_sampling_on_photon_emit;
    
public:
    RENDERING_RECTANGLE_LIGHT(const VECTOR<T,3>& position_input,const VECTOR<T,3>& color_input,
        const T brightness_input,const VECTOR<T,3>& u_direction_input,const VECTOR<T,3>& v_direction_input,
        const int u_samples_input,const int v_samples_input,RENDER_WORLD<T>& world_input,
        const bool supports_global_photons,const bool supports_caustic_photons,
        const bool supports_volume_photons,const bool use_stratified_sampling_on_photon_emit);
    virtual ~RENDERING_RECTANGLE_LIGHT();

    void Sample_Points(const VECTOR<T,3>& surface_position,const VECTOR<T,3>& surface_normal,
        ARRAY<RAY<VECTOR<T,3> > >& sample_array)const override;
    VECTOR<T,3> Emitted_Light(const RENDERING_RAY<T>& ray) const override;
    int Emit_Photons(RENDERING_RAY<T>& parent_ray,PHOTON_MAP<T>& photon_map,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type) const  override;
//#####################################################################
};
}
#endif
