//#####################################################################
// Copyright 2002, Ronald Fedkiw, Eilene Hao, Sergey Koltakov, Neil Molino, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_BEZIER_SPLINE_2D
//##################################################################### 
#ifndef __OPENGL_BEZIER_SPLINE_2D__
#define __OPENGL_BEZIER_SPLINE_2D__

#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T,int d>
class OPENGL_BEZIER_SPLINE_2D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,d+1> CV;
public:
    const BEZIER_SPLINE<TV,d>& curve;
    OPENGL_COLOR color;
    OPENGL_COLOR vertex_color,vertex_position_color,velocity_color;
    bool draw_vertices,draw_velocities;
    T velocity_scale;
    int selected_vertex;
    int selected_segment;

    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    OPENGL_BEZIER_SPLINE_2D(STREAM_TYPE stream_type,const BEZIER_SPLINE<TV,d>& curve_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan());

    void Display() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    virtual void Clear_Selection() override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;
    RANGE<VECTOR<T,3> > Bounding_Box() const override;
    void Print_Selection_Info(std::ostream& output_stream) const override;
    void Print_Selection_Info(std::ostream& output_stream,MATRIX<T,3>* transform) const;
    void Draw_Highlighted_Spline(int id) const;

    TV Evaluate(int id,T t) const;

private:
    void Draw_Vertices_For_Selection() const;
    void Draw_Curves_For_Selection() const;
};

}
#endif
