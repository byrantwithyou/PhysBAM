//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_POINTS_2D
//#####################################################################
#ifndef __OPENGL_POINTS_2D__
#define __OPENGL_POINTS_2D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Vectors/VECTOR_2D.h>
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
    T point_size;
    bool draw_point_numbers,draw_radii; 
    ARRAY<OPENGL_COLOR>* point_colors;
    ARRAY<int>* point_ids;
    ARRAY<T>* point_radii;

    int selected_index;
    int selected_id;
    OPENGL_COLOR selected_old_color;

    OPENGL_POINTS_2D(T_ARRAY& points_input,const OPENGL_COLOR &color_input = OPENGL_COLOR::White(),float point_size = 5);
    virtual ~OPENGL_POINTS_2D();

    void Set_Points_From_Particles(const GEOMETRY_PARTICLES<TV>& particles,bool keep_colors=true);

    bool Use_Bounding_Box() const override {return points.Size()>0;}
    virtual int Particle_Index(const int index) const {return index;}
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;
    void Display() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    void Print_Selection_Info(std::ostream &output_stream) const override;

    void Store_Point_Colors(bool store_point_colors = true);
    void Store_Point_Ids(bool store_ids=true);
    void Store_Point_Radii(bool store_point_radii = true);

    void Set_Point_Color(int index,const OPENGL_COLOR& point_color);
    void Set_Point_Colors(const ARRAY<int>& indices,const OPENGL_COLOR& point_color);
    void Reset_Point_Colors();

    void Select_Point(int index);
    void Select_Points(const ARRAY<int> &indices);
};

}

#endif
