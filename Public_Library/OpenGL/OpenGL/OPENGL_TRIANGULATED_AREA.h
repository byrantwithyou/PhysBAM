//#####################################################################
// Copyright 2001-2004, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TRIANGULATED_AREA
//##################################################################### 
#ifndef __OPENGL_TRIANGULATED_AREA__
#define __OPENGL_TRIANGULATED_AREA__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T> class TRIANGULATED_AREA;

template<class T>
class OPENGL_TRIANGULATED_AREA:public OPENGL_OBJECT<T>
{
public:
    typedef VECTOR<T,2> TV;
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Point;
    using OPENGL_OBJECT<T>::World_Space_Box;
    TRIANGULATED_AREA<T>& triangulated_area;
    OPENGL_COLOR vertex_color,segment_color,triangle_color,triangle_inverted_color,velocity_color;
    int selected_vertex;
    int selected_segment;
    int selected_triangle;
    ARRAY<OPENGL_COLOR>* color_map;
    bool draw_vertices,draw_velocities;
    T velocity_scale;

    OPENGL_TRIANGULATED_AREA(STREAM_TYPE stream_type,TRIANGULATED_AREA<T>& triangulated_area_input,const bool draw_vertices_input=false,
                             const OPENGL_COLOR& vertex_color_input=OPENGL_COLOR::Red(),
                             const OPENGL_COLOR& segment_color_input=OPENGL_COLOR::Green(),
                             const OPENGL_COLOR& triangle_color_input=OPENGL_COLOR::Blue(),
                             const OPENGL_COLOR& triangle_inverted_color_input=OPENGL_COLOR::Violet(0.75),
                             ARRAY<OPENGL_COLOR>* color_map_input=0);

    void Display() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    void Print_Selection_Info(std::ostream& output_stream) const override;
    void Print_Selection_Info(std::ostream& output_stream,MATRIX<T,3>* transform) const;
    virtual void Set_Color_Map(ARRAY<OPENGL_COLOR>* color_map_input){color_map=color_map_input;}

    void Get_Vertex_Selection(int index);
    void Get_Segment_Selection(int index);
    void Get_Triangle_Selection(int index);

protected:
    void Draw_Vertices() const;
    void Draw_Segments() const;
    void Draw_Triangles(const bool use_color_map=true) const;
    void Draw_Vertices_For_Selection() const;
    void Draw_Segments_For_Selection() const;
    void Draw_Triangles_For_Selection() const;
};
}
#endif
