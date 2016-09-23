//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_UNIFORM_VOXELS
//#####################################################################
#include <Core/Log/LOG.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Computations/SMOOTH_UNIFORM.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_UNIFORM_VOXELS.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_VOXEL_SHADER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_UNIFORM_VOXELS<T>::
RENDERING_UNIFORM_VOXELS(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,3> >& data_input,const T volumetric_step)
    :grid(grid_input),coarse_grid(grid_input),volumetric_step(volumetric_step),number_of_smoothing_steps(0)
{
    box=grid_input.domain;
    voxel_light_interpolation=&default_voxel_light_interpolation;voxel_source_interpolation=&default_voxel_source_interpolation;
    data.Append(&data_input);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_UNIFORM_VOXELS<T>::
RENDERING_UNIFORM_VOXELS(GRID<TV>& grid_input,GRID<TV>& coarse_grid_input,ARRAY<T,VECTOR<int,3> >& data_input,
    const T volumetric_step)
    :grid(grid_input),coarse_grid(coarse_grid_input),volumetric_step(volumetric_step),number_of_smoothing_steps(0)
{
    box=grid_input.domain;
    voxel_light_interpolation=&default_voxel_light_interpolation;voxel_source_interpolation=&default_voxel_source_interpolation;
    data.Append(&data_input);
}
//#####################################################################
// Function Print_Parameters
//#####################################################################
template<class T> void RENDERING_UNIFORM_VOXELS<T>::
Print_Parameters()
{
    for(int i=0;i<data.Size();i++){
        LOG::cout<<"data_scale("<<i<<")="<<data_scale(i)<<" data_offset("<<i<<")="<<data_offset(i)<<std::endl;
        LOG::cout<<"data_clamp_low_value("<<i<<")="<<data_clamp_low_value(i)<<" data_lowest_value("<<i<<")="<<data_lowest_value(i)<<std::endl;
        LOG::cout<<"data_clamp_high_value("<<i<<")="<<data_clamp_high_value(i)<<" data_highest_value("<<i<<")="<<data_highest_value(i)<<std::endl;
        LOG::cout<<std::endl;}
}
//#####################################################################
// Function Volumetric_Integration_Step
//#####################################################################
template<class T> T RENDERING_UNIFORM_VOXELS<T>::
Volumetric_Integration_Step(const RAY<VECTOR<T,3> > &ray,const T xi) const
{
    return xi*volumetric_step;
}
//#####################################################################
// Function Source_Term
//#####################################################################
template<class T> T RENDERING_UNIFORM_VOXELS<T>::
Source_Term(const int source_term_index,const VECTOR<T,3>& location) const
{
    T value=voxel_source_interpolation->Clamped_To_Array(grid,*data(source_term_index),location);
    value=data_scale(source_term_index)*(value+data_offset(source_term_index));
    if(data_clamp_low_value(source_term_index)) value=max(value,data_lowest_value(source_term_index));
    if(data_clamp_high_value(source_term_index)) value=min(value,data_highest_value(source_term_index));
    return value;
}
//#####################################################################
// Function Get_Node_Locations
//#####################################################################
template<class T> void RENDERING_UNIFORM_VOXELS<T>::
Get_Node_Locations(ARRAY<VECTOR<T,3> >& locations)
{
    locations.Resize(coarse_grid.counts.Product());map_from_accessor_index_to_my_index.Resize(locations.m);int index=0;
    for(int i=0;i<coarse_grid.counts.x;i++)
        for(int j=0;j<coarse_grid.counts.y;j++)
            for(int ij=0;ij<coarse_grid.counts.z;ij++){
                map_from_accessor_index_to_my_index(index)=VECTOR<int,3>(i,j,ij);
                locations(index)=coarse_grid.X(TV_INT(i,j,ij));index++;}
}
//#####################################################################
// Function Use_Precomputed_Light_Data
//#####################################################################
template<class T> bool RENDERING_UNIFORM_VOXELS<T>::
Use_Precomputed_Light_Data(const VECTOR<T,3>& location,const int light_index) const
{
    if(!precompute_single_scattering)return false;
    VECTOR<int,3> index=INTERPOLATION_UNIFORM<TV,VECTOR<T,3> >::Clamped_Index_End_Minus_One(coarse_grid,*precomputed_light(light_index),location);
    int i=index.x,j=index.y,ij=index.z;
    if(coarse_grid.Outside(location))return false;
    bool i_j_ij=(*precomputed_light_valid(light_index))(i,j,ij),i_j_ij1=(*precomputed_light_valid(light_index))(i,j,ij+1),i_j1_ij=(*precomputed_light_valid(light_index))(i,j+1,ij),
        i_j1_ij1=(*precomputed_light_valid(light_index))(i,j+1,ij+1),i1_j_ij=(*precomputed_light_valid(light_index))(i+1,j,ij),i1_j_ij1=(*precomputed_light_valid(light_index))(i+1,j,ij+1),
        i1_j1_ij=(*precomputed_light_valid(light_index))(i+1,j+1,ij),i1_j1_ij1=(*precomputed_light_valid(light_index))(i+1,j+1,ij+1);
    return i_j_ij&&i_j_ij1&&i_j1_ij&&i_j1_ij1&&i1_j_ij&&i1_j_ij1&&i1_j1_ij&&i1_j1_ij1;
}
//#####################################################################
// Function Set_Precomputed_Light_Data
//#####################################################################
template<class T> void RENDERING_UNIFORM_VOXELS<T>::
Set_Precomputed_Light_Data(const int location_index,const int light_index,const VECTOR<T,3>& light_value)
{
    (*precomputed_light(light_index))(map_from_accessor_index_to_my_index(location_index))=light_value;
}
//#####################################################################
// Function Set_Precomputed_Light_Valid
//#####################################################################
template<class T> void RENDERING_UNIFORM_VOXELS<T>::
Set_Precomputed_Light_Valid(const int location_index,const int light_index,const bool value)
{
    VECTOR<int,3> index=map_from_accessor_index_to_my_index(location_index);
    (*precomputed_light_valid(light_index))(index)=value;
}
//#####################################################################
// Function Precomputed_Light_Data
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_UNIFORM_VOXELS<T>::
Precomputed_Light_Data(const VECTOR<T,3>& location,const int light) const
{
    return voxel_light_interpolation->Clamped_To_Array(coarse_grid,*precomputed_light(light),location);
}
//#####################################################################
// Function Set_Custom_Source_Interpolation
//#####################################################################
template<class T> void RENDERING_UNIFORM_VOXELS<T>::
Set_Custom_Source_Interpolation(INTERPOLATION_UNIFORM<TV,T>* interpolation)
{
    voxel_source_interpolation=interpolation;
}
//#####################################################################
// Function Set_Custom_Light_Interpolation
//#####################################################################
template<class T> void RENDERING_UNIFORM_VOXELS<T>::
Set_Custom_Light_Interpolation(INTERPOLATION_UNIFORM<TV,VECTOR<T,3> >* interpolation)
{
    voxel_light_interpolation=interpolation;
}
//#####################################################################
// Function Prepare_For_Precomputation
//#####################################################################
template<class T> void RENDERING_UNIFORM_VOXELS<T>::
Prepare_For_Precomputation(RENDER_WORLD<T>& world)
{
    precomputed_light.Resize(world.Lights().m);precomputed_light_valid.Resize(world.Lights().m);
    for(int i=0;i<precomputed_light.m;i++)
        precomputed_light(i)=new ARRAY<VECTOR<T,3> ,VECTOR<int,3> >(0,coarse_grid.counts.x,0,coarse_grid.counts.y,0,coarse_grid.counts.z);
    for(int i=0;i<precomputed_light.m;i++){
        precomputed_light_valid(i)=new ARRAY<bool,VECTOR<int,3> >(0,coarse_grid.counts.x,0,coarse_grid.counts.y,0,coarse_grid.counts.z);
        precomputed_light_valid(i)->Fill(false);}
}
//#####################################################################
// Function Postprocess_Light_Field
//#####################################################################
template<class T> void RENDERING_UNIFORM_VOXELS<T>::
Postprocess_Light_Field()
{
    if(number_of_smoothing_steps) for(int light=0;light<precomputed_light.m;light++){
            LOG::cout<<"Smoothing light "<<light<<" "<<number_of_smoothing_steps<<" steps"<<std::endl;
            SMOOTH::Smooth<T,3>(*precomputed_light(light),number_of_smoothing_steps,0);}
}
//#####################################################################
template class RENDERING_UNIFORM_VOXELS<float>;
template class RENDERING_UNIFORM_VOXELS<double>;
}
