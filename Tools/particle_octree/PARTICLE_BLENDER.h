//#####################################################################
// Copyright 2002-2005, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_BLENDER
//##################################################################### 
#ifndef __PARTICLE_BLENDER__
#define __PARTICLE_BLENDER__

#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
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
    
    BOX<TV> Get_Bounding_Box(T radius,const VECTOR_3D<T>& position) const
    {
        T radius_of_influence=R*radius;
        VECTOR_3D<T> radius_vector(radius_of_influence,radius_of_influence,radius_of_influence);
        return BOX<TV>(position-radius_vector,position+radius_vector);
    }

    BOX<TV> Get_Bounding_Box(T radius_x,T radius_yz,const VECTOR_3D<T>& position,const VECTOR_3D<T>& velocity) const
    {
        if(velocity.Magnitude_Squared()==0) return Get_Bounding_Box(radius_x,position);
        MATRIX_4X4<T> transform=MATRIX_4X4<T>::Rotation_Matrix(VECTOR_3D<T>(1,0,0),velocity);
        T radius_of_influence_x=R*radius_x,radius_of_influence_yz=R*radius_yz;
        ORIENTED_BOX<TV> oriented_bounding_box(BOX<TV>(-radius_of_influence_x,radius_of_influence_x,-radius_of_influence_yz,radius_of_influence_yz,-radius_of_influence_yz,radius_of_influence_yz),transform);
        oriented_bounding_box.corner+=position;
        BOX<TV> aligned_box=oriented_bounding_box.Axis_Aligned_Bounding_Box();
        return aligned_box;
    }

    T Get_Distance(T radius,const VECTOR_3D<T>& position,VECTOR_3D<T>& location) const
    {
        return (location-position).Magnitude()/radius;
    }

    T Get_Distance(T radius_x,T radius_yz,const VECTOR_3D<T>& position,const VECTOR_3D<T>& velocity,VECTOR_3D<T>& location) const
    {
        if(velocity.Magnitude_Squared()==0) return Get_Distance(radius_x,position,location);
        MATRIX_4X4<T> rotation=MATRIX_4X4<T>::Rotation_Matrix(velocity,VECTOR_3D<T>(1,0,0));
        return Get_Distance(radius_x,radius_yz,position,velocity,location,rotation);
    }

    T Get_Distance(T radius_x,T radius_yz,const VECTOR_3D<T>& position,const VECTOR_3D<T>& velocity,VECTOR_3D<T>& location,const MATRIX_4X4<T>& rotation) const
    {
        VECTOR_3D<T> transformed_location=rotation*(location-position);
        return (VECTOR_3D<T>(transformed_location.x/radius_x,transformed_location.y/radius_yz,transformed_location.z/radius_yz).Magnitude());
    }
//#####################################################################
//#####################################################################
};
}
#endif
