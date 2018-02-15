//#####################################################################
// Copyright 2002-2005, Robert Bridson, Eilene Hao, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_AXES
//##################################################################### 
#ifndef __OPENGL_AXES__
#define __OPENGL_AXES__

#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/FRAME.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
namespace PhysBAM{

template<class T>
class OPENGL_AXES:public OPENGL_OBJECT<T>
{
public:
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    using OPENGL_OBJECT<T>::frame;using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;
    using OPENGL_OBJECT<T>::World_Space_Box;
    RANGE<TV> box; // extents of axes with respect to local frame
    bool draw_box; // whether to draw a bounding box or axis vectors
    VECTOR<bool,TV::SPIN::m> draw_grid;
    int number_grid_spaces;

    OPENGL_AXES(const FRAME<TV>& frame_input=FRAME<TV>(),const RANGE<TV>& box_input=RANGE<TV>::Unit_Box());

    virtual ~OPENGL_AXES();

    void Scale(const T scale)
    {box*=scale;}

//#####################################################################
    void Display() const override;
    RANGE<TV> Bounding_Box() const override;
    bool Is_Transparent() const override {return false;}
//#####################################################################
};
}
#endif
