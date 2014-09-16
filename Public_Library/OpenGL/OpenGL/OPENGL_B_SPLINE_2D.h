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
class OPENGL_B_SPLINE_2D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,d+1> CV;
public:
    const B_SPLINE<TV,d>& curve;
    mutable BEZIER_SPLINE<TV,d> bezier_version;
    OPENGL_BEZIER_SPLINE_2D<T,3> bezier_opengl;
    OPENGL_COLOR color; // Currently the rest of these aren't used! So attempts to change them do nothing.
    OPENGL_COLOR vertex_color,vertex_position_color,velocity_color;
    bool draw_vertices,draw_velocities;
    T velocity_scale;
private:
    OPENGL_SELECTION<T>* current_selection;

public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    OPENGL_B_SPLINE_2D(const B_SPLINE<TV,d>& curve_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan());

    void Display() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* selection,MATRIX<T,3>* transform) const;
    void Draw_Highlighted_Spline(int id) const;

    OPENGL_SELECTION<T>* Get_Vertex_Selection(int index);
    OPENGL_SELECTION<T>* Get_Segment_Selection(int index);
    TV Evaluate(int id,T t) const;

//private:
//    void Draw_Vertices_For_Selection() const;
//    void Draw_Curves_For_Selection() const;
};

template<class T,int d>
class OPENGL_SELECTION_B_SPLINE_VERTEX_2D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_B_SPLINE_VERTEX_2D(OPENGL_OBJECT<T>* object, int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::B_SPLINE_VERTEX_2D, object), index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T,int d>
class OPENGL_SELECTION_B_SPLINE_SEGMENT_2D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_B_SPLINE_SEGMENT_2D(OPENGL_OBJECT<T>* object, int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::B_SPLINE_SEGMENT_2D, object), index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}
#endif
