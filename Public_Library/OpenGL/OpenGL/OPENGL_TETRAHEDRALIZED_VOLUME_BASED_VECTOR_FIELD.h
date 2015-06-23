//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD
//#####################################################################
#ifndef __OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD__
#define __OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;
    
template<class T_input>
class OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD:public OPENGL_VECTOR_FIELD_3D<T_input>
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    using OPENGL_VECTOR_FIELD_3D<T>::vector_field;using OPENGL_VECTOR_FIELD_3D<T>::vector_locations;
    using OPENGL_VECTOR_FIELD_3D<T>::size;using OPENGL_OBJECT<T>::World_Space_Box;

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume;
    ARRAY<TV>& V;

    OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD(STREAM_TYPE stream_type,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,ARRAY<TV>& V);
    virtual ~OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD();

    void Update();  // Call when tetrahedralized volume/V change

    virtual RANGE<TV> Bounding_Box() const override;
};
}
#endif
