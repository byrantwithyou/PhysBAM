//#####################################################################
// Copyright 2002, 2003, Ronald Fedkiw, Eran Guendelman, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_VECTOR_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_VECTOR_FIELD_3D__
#define __OPENGL_VECTOR_FIELD_3D__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T>
class OPENGL_VECTOR_FIELD_3D:public OPENGL_OBJECT<T>
{
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Point;
    ARRAY<VECTOR<T,3> >& vector_field;
    ARRAY<VECTOR<T,3> >& vector_locations;
    OPENGL_COLOR vector_color;
    double size;
    bool draw_arrowhead;
    bool draw_value;
    bool draw_basepoint;
    bool draw_fancy_arrow;
    mutable GLUquadric* vector_hat;

    OPENGL_VECTOR_FIELD_3D(ARRAY<VECTOR<T,3> >& field,ARRAY<VECTOR<T,3> >& locations, 
        const OPENGL_COLOR& color=OPENGL_COLOR::White(),double size=0.025, 
        bool draw_arrowhead=false,bool draw_value=false,bool draw_basepoint=false,bool draw_fancy_arrow=false);

    virtual ~OPENGL_VECTOR_FIELD_3D();

//#####################################################################
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Display() const PHYSBAM_OVERRIDE;
    void Scale_Vector_Size(const T scale);
    void Toggle_Arrowhead_Mode();
//#####################################################################
};
}
#endif
