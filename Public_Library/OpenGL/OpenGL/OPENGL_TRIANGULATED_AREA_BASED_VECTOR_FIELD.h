//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD
//#####################################################################
#ifndef __OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD__
#define __OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
namespace PhysBAM{

template<class T> class TRIANGULATED_AREA;
    
template<class T_input>
class OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD:public OPENGL_VECTOR_FIELD_2D<ARRAY<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_OBJECT<T>::World_Space_Point;using OPENGL_OBJECT<T>::World_Space_Box;
    typedef OPENGL_VECTOR_FIELD_2D<ARRAY<TV> > BASE;
    using BASE::size;

    TRIANGULATED_AREA<T>& triangulated_area;
    ARRAY<TV>& V;
    ARRAY<TV> vector_field,vector_locations;

    OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD(TRIANGULATED_AREA<T>& triangulated_area,ARRAY<TV>& V);
    virtual ~OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD();

    void Update();  // Call when triangulated area/V change

    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif
