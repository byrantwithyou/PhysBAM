//#####################################################################
// Class OPENGL_B_SPLINE_2D
//##################################################################### 
#ifndef __OPENGL_B_SPLINE_2D__
#define __OPENGL_B_SPLINE_2D__

#include <Geometry/Topology_Based_Geometry/B_SPLINE.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <OpenGL/OpenGL/OPENGL_BEZIER_SPLINE_2D.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T,int d>
class OPENGL_B_SPLINE_2D:public OPENGL_BEZIER_SPLINE_2D<T,d>
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,d+1> CV;
public:
    typedef OPENGL_BEZIER_SPLINE_2D<T,d> BASE;
    mutable BEZIER_SPLINE<TV,d> bezier_version;
    const B_SPLINE<TV,d>& curve;

public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    OPENGL_B_SPLINE_2D(const B_SPLINE<TV,d>& curve_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan());

    void Display() const override;
};

}
#endif
