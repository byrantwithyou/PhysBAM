//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_UNIFORM_VOXELS
//#####################################################################
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Grids_Uniform_Computations/SMOOTH_UNIFORM.h>
#include <Tools/Log/LOG.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_UNIFORM_VOXELS.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_VOXEL_SHADER.h>
using namespace PhysBAM;
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
namespace PhysBAM{
template class RENDERING_UNIFORM_VOXELS<float>;
template class RENDERING_UNIFORM_VOXELS<double>;
}
