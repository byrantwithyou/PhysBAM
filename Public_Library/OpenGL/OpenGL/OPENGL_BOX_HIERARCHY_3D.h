//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_BOX_HIERARCHY_3D
//##################################################################### 
#ifndef __OPENGL_BOX_HIERARCHY_3D__
#define __OPENGL_BOX_HIERARCHY_3D__

#include <Tools/Math_Tools/RANGE.h>
#include <Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T>
class OPENGL_BOX_HIERARCHY_3D:public OPENGL_OBJECT
{
public:
    BOX_HIERARCHY<VECTOR<T,3> >* hierarchy;
    OPENGL_COLOR color;
    int min_height,max_height;

    OPENGL_BOX_HIERARCHY_3D(BOX_HIERARCHY<VECTOR<T,3> > *hierarchy_in,const OPENGL_COLOR &color_input = OPENGL_COLOR::White()) 
        :hierarchy(hierarchy_in),color(color_input),min_height(1),max_height(1)
    {}

    void Display() const PHYSBAM_OVERRIDE;
    void Display_Helper(const int cell,const int height) const;
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Increment_Height();
    void Decrement_Height();

//##################################################################### 
};
}
#endif
