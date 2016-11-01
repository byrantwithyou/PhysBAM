//#####################################################################
// Copyright 2002, Ronald Fedkiw, Eilene Hao, Sergey Koltakov, Neil Molino, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SEGMENTED_CURVE_2D
//##################################################################### 
#ifndef __OPENGL_SEGMENTED_CURVE_2D__
#define __OPENGL_SEGMENTED_CURVE_2D__

#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T>
class OPENGL_SEGMENTED_CURVE_2D:public OPENGL_OBJECT<T>
{
public:
    typedef VECTOR<T,2> TV;
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    const SEGMENTED_CURVE_2D<T>& curve;
    OPENGL_COLOR color;
    OPENGL_COLOR vertex_color,vertex_position_color,velocity_color;
    bool draw_vertices,draw_velocities;
    T velocity_scale;
    int selected_vertex;
    int selected_segment;

public:
    OPENGL_SEGMENTED_CURVE_2D(STREAM_TYPE stream_type,const SEGMENTED_CURVE_2D<T>& curve_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan());

    void Display() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    void Print_Selection_Info(std::ostream &output_stream) const override;
    void Print_Selection_Info(std::ostream& output_stream,MATRIX<T,3>* transform) const;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;

private:
    void Draw_Vertices_For_Selection() const;
    void Draw_Segments_For_Selection() const;
};

}
#endif
