//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/RANGE.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>::
OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD(TRIANGULATED_AREA<T>& triangulated_area,ARRAY<TV>& V)
    :OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >(vector_field,vector_locations),triangulated_area(triangulated_area),V(V)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>::
~OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD()
{}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>::
Update()
{
    V.Resize(triangulated_area.particles.Size());
    vector_field.Resize(V.m);
    vector_locations.Resize(V.m);
    for(int i=0;i<triangulated_area.particles.Size();i++){
        vector_field(i)=V(i);vector_locations(i)=triangulated_area.particles.X(i);}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>::
Bounding_Box() const
{
    if(triangulated_area.bounding_box) return World_Space_Box(*triangulated_area.bounding_box);
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<float>;
template class OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<double>;
}
