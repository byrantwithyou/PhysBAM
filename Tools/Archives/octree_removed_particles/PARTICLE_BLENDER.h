//#####################################################################
// Copyright 2002-2005, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_BLENDER
//##################################################################### 
#ifndef __PARTICLE_BLENDER__
#define __PARTICLE_BLENDER__

#include <PhysBAM_Tools/Polynomials/CUBIC.h>

namespace PhysBAM{

template<class T>
class PARTICLE_BLENDER
{
    T blending_parameter;
    CUBIC<T> kernel;
    T small_number;
    T R;

public:
    PARTICLE_BLENDER(T blending_parameter_input)
        :blending_parameter(blending_parameter_input),kernel(2,-3,0,1),small_number(1e-8)
    {
        blending_parameter=(blending_parameter>=1)?1-small_number:blending_parameter;
        blending_parameter=(blending_parameter<=0)?small_number:blending_parameter;
        kernel.c0-=blending_parameter;
        kernel.Compute_Roots_Noniterative_In_Interval(0,1);R=1/kernel.root1;
        std::cout<<"Roots:"<<kernel.roots<<std::endl;
        std::cout<<"root1:"<<kernel.root1<<std::endl;
        if(kernel.roots<1){std::cerr<<"Error: kernel.roots=="<<kernel.roots<<std::endl;exit(1);}
        kernel.c0=1;
    }

    T Get_R() const {return R;}

    T C(const T r) const
    {
        if(r>R) return 0;        
        return kernel(r/R);
    }

    BOX_3D<T> Get_Bounding_Box(T radius_x,T radius_yz,const VECTOR<T,3>& position,const VECTOR<T,3>& major_axis) const
    {
        T radius_of_influence_x=R*radius_x,radius_of_influence_yz=R*radius_yz;
        VECTOR<T,3> orthogonal_vector=major_axis.Orthogonal_Vector().Robust_Normalized();
        VECTOR<T,3> scaled_x_axis=radius_of_influence_x*major_axis;
        VECTOR<T,3> scaled_y_axis=orthogonal_vector*radius_of_influence_yz;
        VECTOR<T,3> scaled_z_axis=VECTOR<T,3>::Cross_Product(major_axis,orthogonal_vector).Normalized()*radius_of_influence_yz;
        ORIENTED_BOX_3D<T> oriented_bounding_box(position-scaled_x_axis-scaled_y_axis-scaled_z_axis,(T)2*scaled_x_axis,(T)2*scaled_y_axis,(T)2*scaled_z_axis);
        return oriented_bounding_box.Axis_Aligned_Bounding_Box();
    }

    T Get_Distance(double one_over_radius_x_squared,double one_over_radius_yz_squared,const VECTOR<T,3>& position,const VECTOR<T,3>& major_axis,VECTOR<T,3>& location) const
    {
        VECTOR<T,3> dX=location-position;T dot=VECTOR<T,3>::Dot_Product(dX,major_axis);
        return sqrt(one_over_radius_yz_squared*dX.Magnitude_Squared() + (one_over_radius_x_squared-one_over_radius_yz_squared)*sqr(dot));
    }
//#####################################################################
//#####################################################################
};
}
#endif
