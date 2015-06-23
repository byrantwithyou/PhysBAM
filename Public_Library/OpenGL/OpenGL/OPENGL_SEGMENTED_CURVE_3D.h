//#####################################################################
// Copyright 2002-2007, Kevin Der, Ronald Fedkiw, Eilene Hao, Sergey Koltakov, Michael Lentine, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SEGMENTED_CURVE_3D
//##################################################################### 
#ifndef __OPENGL_SEGMENTED_CURVE_3D__
#define __OPENGL_SEGMENTED_CURVE_3D__

#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T> class OPENGL_SELECTION;

template<class T>
class OPENGL_SEGMENTED_CURVE_3D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    const SEGMENTED_CURVE<TV>& curve;
    const OPENGL_SEGMENTED_CURVE_3D<T>* parent_curve;
    mutable ARRAY<int> segment_nodes;
    mutable HASHTABLE<int,TV> vertex_normals;
    OPENGL_COLOR color;
    OPENGL_COLOR vertex_color,vertex_position_color;
    bool draw_vertices,use_solid_color,hide_unselected;

    OPENGL_SEGMENTED_CURVE_3D(STREAM_TYPE stream_type,const SEGMENTED_CURVE<TV>& curve_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan())
        :OPENGL_OBJECT<T>(stream_type),curve(curve_input),parent_curve(0),color(color_input),vertex_color(OPENGL_COLOR::Green(0.9)),
        vertex_position_color(OPENGL_COLOR::Magenta()),draw_vertices(false),use_solid_color(true),smooth_normals(false),current_selection(0)
    {}

    void Display() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer,int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) override;
    void Clear_Highlight() override;
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const override;

    OPENGL_SELECTION<T>* Get_Vertex_Selection(int index);
    OPENGL_SELECTION<T>* Get_Segment_Selection(int index);
    OPENGL_SELECTION<T>* Get_Curve_Selection(int index);

    void Turn_Smooth_Shading_On() override;
    void Turn_Smooth_Shading_Off() override;

    void Initialize_Vertex_Normals() const;
    ARRAY<int> Get_Selected_Edges() const;

private:
    void Draw_Vertices_For_Selection() const;
    void Draw_Segments_For_Selection() const;

    bool smooth_normals;

    OPENGL_SELECTION<T>* current_selection;
};

template<class T>
class OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D(OPENGL_OBJECT<T>* object,int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::SEGMENTED_CURVE_VERTEX_3D,object),index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

template<class T>
class OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D(OPENGL_OBJECT<T>* object,int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_3D,object),index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

template<class T>
class OPENGL_SELECTION_SEGMENTED_CURVE_3D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    OPENGL_SELECTION_SEGMENTED_CURVE_3D(OPENGL_OBJECT<T>* object,int index=0) 
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::SEGMENTED_CURVE_3D,object),index(index) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

}
#endif
