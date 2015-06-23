//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_POINT_SIMPLICES_1D
//##################################################################### 
#ifndef __OPENGL_POINT_SIMPLICES_1D__
#define __OPENGL_POINT_SIMPLICES_1D__

#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T> class OPENGL_SELECTION;

template<class T>
class OPENGL_POINT_SIMPLICES_1D:public OPENGL_OBJECT<T>
{
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    const POINT_SIMPLICES_1D<T>& simplices;
    OPENGL_COLOR color;
    OPENGL_COLOR vertex_color,segment_color,vertex_position_color,velocity_color;
    bool draw_vertices;

    OPENGL_POINT_SIMPLICES_1D(STREAM_TYPE stream_type,const POINT_SIMPLICES_1D<T>& simplices_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan());

//#####################################################################
    void Display() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
//#####################################################################
};
}
#endif
