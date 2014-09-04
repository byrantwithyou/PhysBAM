//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_FACE_SCALAR_FIELD_2D
//#####################################################################
#ifndef __OPENGL_FACE_SCALAR_FIELD_2D__
#define __OPENGL_FACE_SCALAR_FIELD_2D__

#include <Tools/Grids_Uniform/GRID.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_POINTS_2D.h>

namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_FACE_SCALAR_FIELD_2D:public OPENGL_OBJECT<T>
{
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;GRID<TV> grid;
    ARRAY<T2,FACE_INDEX<2> > &face_values;
    ARRAY_VIEW<T2,VECTOR<int,2> > &x_face_values,&y_face_values;
    OPENGL_COLOR_MAP<T2> *color_map;
private:
    OPENGL_POINTS_2D<T> opengl_points;

//#####################################################################
public:
    OPENGL_FACE_SCALAR_FIELD_2D(const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<2> > &face_values_input,OPENGL_COLOR_MAP<T2> *color_map_input);
    virtual ~OPENGL_FACE_SCALAR_FIELD_2D();
    void Display() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual void Update();  // Call when values or other attributes have changed
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
