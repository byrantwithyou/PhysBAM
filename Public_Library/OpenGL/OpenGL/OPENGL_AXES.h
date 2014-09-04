//#####################################################################
// Copyright 2002-2005, Robert Bridson, Eilene Hao, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_AXES
//##################################################################### 
#ifndef __OPENGL_AXES__
#define __OPENGL_AXES__

#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/FRAME.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
namespace PhysBAM{

template<class T>
class OPENGL_AXES:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_OBJECT<T>::frame;using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;
    using OPENGL_OBJECT<T>::World_Space_Box;
    RANGE<TV> box; // extents of axes with respect to local frame
    bool draw_box; // whether to draw a bounding box or axis vectors
    bool draw_xz_grid,draw_xy_grid,draw_yz_grid; // whether to draw grids on each plane
    T grid_spacing;

    OPENGL_AXES(const FRAME<TV>& frame_input=FRAME<TV>(),const RANGE<TV>& box_input=RANGE<TV>::Unit_Box(),
        bool draw_box_input=false,bool draw_xz_grid_input=false,bool draw_xy_grid_input=false,bool draw_yz_grid_input=false,T grid_spacing_input=.1);

    virtual ~OPENGL_AXES();

    void Scale(const T scale)
    {box*=scale;}

//#####################################################################
    void Display() const PHYSBAM_OVERRIDE;
    RANGE<TV> Bounding_Box() const PHYSBAM_OVERRIDE;
    bool Is_Transparent() const PHYSBAM_OVERRIDE {return false;}
//#####################################################################
};
}
#endif
