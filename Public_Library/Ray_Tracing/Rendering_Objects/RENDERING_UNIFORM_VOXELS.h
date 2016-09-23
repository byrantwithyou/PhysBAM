//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_UNIFORM_VOXELS
//#####################################################################
#ifndef __RENDERING_UNIFORM_VOXELS__
#define __RENDERING_UNIFORM_VOXELS__

#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
namespace PhysBAM{

template<class T>
class RENDERING_UNIFORM_VOXELS:public RENDERING_VOXELS<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    using RENDERING_VOXELS<T>::box;using RENDERING_VOXELS<T>::small_number;using RENDERING_VOXELS<T>::precompute_single_scattering;

    GRID<TV>& grid;
    GRID<TV>& coarse_grid;
    ARRAY<ARRAY<T,VECTOR<int,3> >*> data; // defined at each grid point
    ARRAY<T> data_scale;
    ARRAY<T> data_offset;
    ARRAY<bool> data_clamp_low_value;
    ARRAY<T> data_lowest_value;
    ARRAY<bool> data_clamp_high_value;
    ARRAY<T> data_highest_value;
    ARRAY<ARRAY<VECTOR<T,3> ,VECTOR<int,3> >*> precomputed_light;
    ARRAY<ARRAY<bool,VECTOR<int,3> >*> precomputed_light_valid;
    ARRAY<VECTOR<int,3> > map_from_accessor_index_to_my_index;
    T volumetric_step;
    INTERPOLATION_UNIFORM<TV,T>* voxel_source_interpolation;
    INTERPOLATION_UNIFORM<TV,VECTOR<T,3> >* voxel_light_interpolation;
    LINEAR_INTERPOLATION_UNIFORM<TV,T> default_voxel_source_interpolation;
    LINEAR_INTERPOLATION_UNIFORM<TV,VECTOR<T,3> > default_voxel_light_interpolation;
    int number_of_smoothing_steps;

    RENDERING_UNIFORM_VOXELS(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,3> >& data_input,const T volumetric_step);
    RENDERING_UNIFORM_VOXELS(GRID<TV>& grid_input,GRID<TV>& coarse_grid_input,ARRAY<T,VECTOR<int,3> >& data_input,
        const T volumetric_step);

//#####################################################################
    void Print_Parameters();
    T Volumetric_Integration_Step(const RAY<VECTOR<T,3> > &ray,const T xi) const override;
    T Source_Term(const int source_term_index,const VECTOR<T,3>& location) const override;
    void Get_Node_Locations(ARRAY<VECTOR<T,3> >& locations) override;
    bool Use_Precomputed_Light_Data(const VECTOR<T,3>& location,const int light_index) const override;
    void Set_Precomputed_Light_Data(const int location_index,const int light_index,const VECTOR<T,3>& light_value) override;
    void Set_Precomputed_Light_Valid(const int location_index,const int light_index,const bool value) override;
    VECTOR<T,3> Precomputed_Light_Data(const VECTOR<T,3>& location,const int light) const override;
    void Set_Custom_Source_Interpolation(INTERPOLATION_UNIFORM<TV,T>* interpolation);
    void Set_Custom_Light_Interpolation(INTERPOLATION_UNIFORM<TV,VECTOR<T,3> >* interpolation);
protected:
    void Prepare_For_Precomputation(RENDER_WORLD<T>& world) override;
    void Postprocess_Light_Field() override;
//#####################################################################
};
}
#endif
