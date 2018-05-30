//#####################################################################
// Copyright 2002-2005, Eran Guendelman, Geoffrey Irving, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REMOVED_PARTICLES_BLENDER_3D
//##################################################################### 
#include <Core/Log/LOG.h>
#include <Geometry/Basic_Geometry/ELLIPSOID.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Dynamics/Level_Sets/REMOVED_PARTICLES_BLENDER_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> REMOVED_PARTICLES_BLENDER_3D<T>::
REMOVED_PARTICLES_BLENDER_3D(T blending_parameter_input)
    :blending_parameter(blending_parameter_input),kernel(2,-3,0,1),small_number((T)1e-8)
{
    blending_parameter=clamp(blending_parameter,small_number,1-small_number);
    kernel.c0-=blending_parameter;
    kernel.Compute_Roots_Noniterative_In_Interval(0,1);R=1/kernel.root[0];
    LOG::cout<<"roots = "<<kernel.roots<<std::endl;
    LOG::cout<<"root[0] = "<<kernel.root[0]<<std::endl;
    if(kernel.roots<1){LOG::cerr<<"Error: kernel.roots=="<<kernel.roots<<std::endl;PHYSBAM_FATAL_ERROR();}
    kernel.c0=1;
}
//#####################################################################
// Function Get_Oriented_Bounding_Box
//#####################################################################
template<class T> ORIENTED_BOX<VECTOR<T,3> > REMOVED_PARTICLES_BLENDER_3D<T>::
Get_Oriented_Bounding_Box(T radius_x,T radius_yz,const TV& position,const TV& major_axis) const
{
    T radius_of_influence_x=R*radius_x,radius_of_influence_yz=R*radius_yz;
    TV orthogonal_vector=major_axis.Unit_Orthogonal_Vector();
    TV scaled_x_axis=radius_of_influence_x*major_axis;
    TV scaled_y_axis=orthogonal_vector*radius_of_influence_yz;
    TV scaled_z_axis=TV::Cross_Product(major_axis,orthogonal_vector).Normalized()*radius_of_influence_yz;
    return ORIENTED_BOX<TV>(position-scaled_x_axis-scaled_y_axis-scaled_z_axis,MATRIX<T,TV::m>(scaled_x_axis,scaled_y_axis,scaled_z_axis)*2);
}
//#####################################################################
// Function Get_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > REMOVED_PARTICLES_BLENDER_3D<T>::
Get_Bounding_Box(T radius_x,T radius_yz,const TV& position,const TV& major_axis) const
{
    return Get_Oriented_Bounding_Box(radius_x,radius_yz,position,major_axis).Axis_Aligned_Bounding_Box();
}
//#####################################################################
// Function Get_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > REMOVED_PARTICLES_BLENDER_3D<T>::
Get_Bounding_Box(const ELLIPSOID<TV>& ellipsoid) const
{
    return ellipsoid.Oriented_Bounding_Box().Scaled_About_Center(R).Axis_Aligned_Bounding_Box();
}
//#####################################################################
namespace PhysBAM{
template class REMOVED_PARTICLES_BLENDER_3D<float>;
template class REMOVED_PARTICLES_BLENDER_3D<double>;
}
