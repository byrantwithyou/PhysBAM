//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_POINTS_2D
//#####################################################################
#ifndef __OPENGL_POINTS_2D__
#define __OPENGL_POINTS_2D__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T,class T_ARRAY=ARRAY<VECTOR<T,2> > >
class OPENGL_POINTS_2D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    T_ARRAY& points;
    OPENGL_COLOR color;
    float point_size;
    bool draw_point_numbers,draw_radii; 
    ARRAY<OPENGL_COLOR>* point_colors;
    ARRAY<int>* point_ids;
    ARRAY<T>* point_radii;

    OPENGL_POINTS_2D(STREAM_TYPE stream_type,T_ARRAY& points_input,const OPENGL_COLOR &color_input = OPENGL_COLOR::White(),float point_size = 5);
    virtual ~OPENGL_POINTS_2D();

    void Set_Points_From_Particles(const GEOMETRY_PARTICLES<TV>& particles,bool keep_colors=true,const bool use_ids=true)
    {points=particles.X;
    const ARRAY_VIEW<int>* id=use_ids?particles.template Get_Array<int>(ATTRIBUTE_ID_ID):0;
    Store_Point_Ids(id!=0);
    if(point_colors && (!keep_colors || point_colors->m!=particles.Size()))
        Store_Point_Colors(false);

    if(const ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR)){
        Store_Point_Colors(true);
        for(int i=0;i<point_colors->m;i++)
            (*point_colors)(i)=OPENGL_COLOR((*color_attribute)(i));}

    if(id) *point_ids=*id;

    if(const ARRAY_VIEW<T>* radius_attribute=particles.template Get_Array<T>(ATTRIBUTE_ID_RADIUS)){
        Store_Point_Radii(true);
        for(int i=0;i<point_radii->m;i++)
            (*point_radii)(i)=(*radius_attribute)(i);}}

    bool Use_Bounding_Box() const override {return points.Size()>0;}
    virtual int Particle_Index(const int index) const {return index;}
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    void Display() const override;

    OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer, int buffer_size) override;
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) override;
    void Clear_Highlight() override;
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const override;

    void Store_Point_Colors(bool store_point_colors = true);
    void Store_Point_Ids(bool store_ids=true);
    void Store_Point_Radii(bool store_point_radii = true);

    void Set_Point_Color(int index,const OPENGL_COLOR &point_color);
    void Set_Point_Colors(const ARRAY<int> &indices,const OPENGL_COLOR &point_color);
    void Reset_Point_Colors();

    void Select_Point(int index);
    void Select_Points(const ARRAY<int> &indices);
    void Clear_Selection();

};

template<class T>
class OPENGL_SELECTION_POINTS_2D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;
    bool has_id;
    int id;

    OPENGL_SELECTION_POINTS_2D(OPENGL_OBJECT<T>* object):OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::POINTS_2D, object) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

}

#endif
