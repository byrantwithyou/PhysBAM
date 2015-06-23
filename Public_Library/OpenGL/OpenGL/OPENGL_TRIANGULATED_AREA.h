//#####################################################################
// Copyright 2001-2004, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TRIANGULATED_AREA
//##################################################################### 
#ifndef __OPENGL_TRIANGULATED_AREA__
#define __OPENGL_TRIANGULATED_AREA__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T> class TRIANGULATED_AREA;
template<class T> class OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX;
template<class T> class OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT;
template<class T> class OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE;

template<class T>
class OPENGL_TRIANGULATED_AREA:public OPENGL_OBJECT<T>
{
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Point;
    TRIANGULATED_AREA<T>& triangulated_area;
    OPENGL_COLOR vertex_color,segment_color,triangle_color,triangle_inverted_color,velocity_color;
    OPENGL_SELECTION<T>* current_selection;
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

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint* buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) override;
    void Clear_Highlight() override;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* selection) const override;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* selection,MATRIX<T,3>* transform) const;
    virtual void Set_Color_Map(ARRAY<OPENGL_COLOR>* color_map_input){color_map=color_map_input;}

    OPENGL_SELECTION<T>* Get_Vertex_Selection(int index);
    OPENGL_SELECTION<T>* Get_Segment_Selection(int index);
    OPENGL_SELECTION<T>* Get_Triangle_Selection(int index);

protected:
    void Draw_Vertices() const;
    void Draw_Segments() const;
    void Draw_Triangles(const bool use_color_map=true) const;
    void Draw_Vertices_For_Selection() const;
    void Draw_Segments_For_Selection() const;
    void Draw_Triangles_For_Selection() const;

    friend class OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T>;
    friend class OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T>;
    friend class OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T>;
};

template<class T>
class OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX(OPENGL_OBJECT<T>* object,int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::TRIANGULATED_AREA_VERTEX,object),index(index){}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

template<class T>
class OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT(OPENGL_OBJECT<T>* object,int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::TRIANGULATED_AREA_SEGMENT,object),index(index){}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

template<class T>
class OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE(OPENGL_OBJECT<T>* object,int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::TRIANGULATED_AREA_TRIANGLE,object),index(index){}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};
}
#endif
