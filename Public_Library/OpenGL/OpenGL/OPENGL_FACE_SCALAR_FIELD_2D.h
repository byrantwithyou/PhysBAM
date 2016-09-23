//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_FACE_SCALAR_FIELD_2D
//#####################################################################
#ifndef __OPENGL_FACE_SCALAR_FIELD_2D__
#define __OPENGL_FACE_SCALAR_FIELD_2D__

#include <Grid_Tools/Grids/GRID.h>
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
    OPENGL_COLOR_MAP<T2> *color_map;
    OPENGL_POINTS_2D<T> opengl_points;

//#####################################################################
    OPENGL_FACE_SCALAR_FIELD_2D(STREAM_TYPE stream_type,const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<2> > &face_values_input,OPENGL_COLOR_MAP<T2> *color_map_input);
    virtual ~OPENGL_FACE_SCALAR_FIELD_2D();
    void Display() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    virtual void Update();  // Call when values or other attributes have changed
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const override;
//#####################################################################
};
}
#endif
