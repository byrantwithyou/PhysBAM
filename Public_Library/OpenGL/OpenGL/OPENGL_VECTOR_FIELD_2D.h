//#####################################################################
// Copyright 2003, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_VECTOR_FIELD_2D
//##################################################################### 
#ifndef __OPENGL_VECTOR_FIELD_2D__
#define __OPENGL_VECTOR_FIELD_2D__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T_ARRAY>
class OPENGL_VECTOR_FIELD_2D:public OPENGL_OBJECT<typename T_ARRAY::SCALAR>
{
    typedef typename T_ARRAY::SCALAR T;
    typedef VECTOR<T,2> TV;
    STATIC_ASSERT(IS_SAME<typename T_ARRAY::ELEMENT,TV>::value);
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Point;
    const T_ARRAY& vector_field;
    const T_ARRAY& vector_locations;
    OPENGL_COLOR vector_color;
    double size;
    bool draw_arrowhead;
    bool draw_value;
    bool draw;

    OPENGL_VECTOR_FIELD_2D(STREAM_TYPE stream_type,T_ARRAY& vector_field,T_ARRAY& vector_locations,const OPENGL_COLOR &color=OPENGL_COLOR::White(),
        double size=0.025,bool draw_arrowhead=true,bool draw_value=false);
    virtual ~OPENGL_VECTOR_FIELD_2D();

//##################################################################### 
    void Display() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Scale_Vector_Size(const T scale);
//##################################################################### 
};
}
#endif
