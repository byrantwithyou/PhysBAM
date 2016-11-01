//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_POINTS_3D
//#####################################################################
#ifndef __OPENGL_POINTS_3D__
#define __OPENGL_POINTS_3D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T,class T_ARRAY=ARRAY<VECTOR<T,3> > >
class OPENGL_POINTS_3D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    T_ARRAY& points;
    OPENGL_COLOR color;
    T point_size;
    bool draw_point_numbers; 
    ARRAY<OPENGL_COLOR>* point_colors;
    ARRAY<int>* point_ids;
    ARRAY<bool>* draw_mask;

    int selected_index;
    int selected_id;
    OPENGL_COLOR selected_old_color;

    OPENGL_POINTS_3D(STREAM_TYPE stream_type,T_ARRAY& points_input,const OPENGL_COLOR& color_input=OPENGL_COLOR::White(),const T point_size=5);
    virtual ~OPENGL_POINTS_3D();

    template<class T_PARTICLES> void Set_Points_From_Particles(const T_PARTICLES& particles,const bool keep_colors=true,const bool use_ids=true,const bool use_draw_mask=false)
    {points.Resize(particles.Size());
    //Store_Point_Ids(use_ids && particles.store_id);
    Store_Point_Colors(false);
    if(use_draw_mask) Store_Draw_Mask(true);
    for(int i=0;i<particles.Size();i++){
        if(draw_mask) (*draw_mask)(i)=true;
        points(i)=particles.X(i);}}

    bool Use_Bounding_Box() const override {return points.Size()>0;}
    virtual int Particle_Index(const int index) const {return index;}
    virtual RANGE<TV> Bounding_Box() const override;
    virtual RANGE<TV> Selection_Bounding_Box() const override;
    void Display() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    void Print_Selection_Info(std::ostream &output_stream) const override;

    void Store_Point_Colors(bool store_point_colors=true);
    void Store_Point_Ids(bool store_ids=true);
    void Store_Draw_Mask(bool store_draw_mask=true);
    
    void Set_Point_Color(int index,const OPENGL_COLOR& point_color);
    void Set_Point_Colors(const ARRAY<int>& indices,const OPENGL_COLOR& point_color);
    void Reset_Point_Colors();

    void Filter_By_Slice_Planes(const PLANE<T>& leftplane,const PLANE<T>& rightplane);
    void Clear_Filter();

    void Select_Point(int index);
    void Select_Points(const ARRAY<int>& indices);
};

}
#endif
