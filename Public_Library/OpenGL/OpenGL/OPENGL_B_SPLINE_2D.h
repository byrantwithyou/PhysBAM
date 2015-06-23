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

template<class T> class OPENGL_SELECTION;
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
    OPENGL_B_SPLINE_2D(STREAM_TYPE stream_type,const B_SPLINE<TV,d>& curve_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan());

    void Display() const override;

    TV Evaluate(int id,T t) const; // What is this line doing here?
};

template<class T,int d>
class OPENGL_SELECTION_B_SPLINE_VERTEX_2D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_B_SPLINE_VERTEX_2D(OPENGL_OBJECT<T>* object, int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::B_SPLINE_VERTEX_2D, object), index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

template<class T,int d>
class OPENGL_SELECTION_B_SPLINE_SEGMENT_2D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_B_SPLINE_SEGMENT_2D(OPENGL_OBJECT<T>* object, int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::B_SPLINE_SEGMENT_2D, object), index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

}
#endif
