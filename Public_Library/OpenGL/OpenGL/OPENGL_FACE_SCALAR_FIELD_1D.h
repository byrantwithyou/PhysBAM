//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_FACE_SCALAR_FIELD_1D
//#####################################################################
#ifndef __OPENGL_FACE_SCALAR_FIELD_1D__
#define __OPENGL_FACE_SCALAR_FIELD_1D__

#include <Tools/Grids_Uniform/GRID.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_FACE_SCALAR_FIELD_1D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,1> TV;
public:
    using OPENGL_OBJECT<T>::World_Space_Box;
    GRID<TV> grid;
    ARRAY<T2,FACE_INDEX<1> > &face_values;
    ARRAY_VIEW<T2,VECTOR<int,1> > &x_face_values;
    OPENGL_COLOR point_color;
    OPENGL_COLOR line_color;
private:
    T scale;
    
public:
    void Scale(const T scale_input)
    {scale=scale*scale_input;}

//#####################################################################
    OPENGL_FACE_SCALAR_FIELD_1D(const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<1> > &face_values_input,OPENGL_COLOR point_color_input,OPENGL_COLOR line_color_input);
    virtual ~OPENGL_FACE_SCALAR_FIELD_1D();
    void Display() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
