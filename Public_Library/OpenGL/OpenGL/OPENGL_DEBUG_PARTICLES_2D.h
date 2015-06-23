//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_DEBUG_PARTICLES_2D
//#####################################################################
#ifndef __OPENGL_DEBUG_PARTICLES_2D__
#define __OPENGL_DEBUG_PARTICLES_2D__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Particles/PARTICLES.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T>
class OPENGL_DEBUG_PARTICLES_2D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<DEBUG_OBJECT<TV> >& debug_objects;
    OPENGL_COLOR default_color;
    OPENGL_COLOR velocity_color;
    bool draw_velocities;
    bool draw_arrows;
    T scale_velocities;

    OPENGL_DEBUG_PARTICLES_2D(STREAM_TYPE stream_type,GEOMETRY_PARTICLES<TV>& particles_input,ARRAY<DEBUG_OBJECT<TV> >& debug_objects_input,const OPENGL_COLOR &color_input = OPENGL_COLOR::White());
    ~OPENGL_DEBUG_PARTICLES_2D();

    bool Use_Bounding_Box() const override;
    RANGE<VECTOR<T,3> > Bounding_Box() const override;
    void Display() const override;

    OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer, int buffer_size) override;
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) override;
    void Clear_Highlight() override;
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const override;

    void Select_Point(int index);
    void Select_Points(const ARRAY<int> &indices);
    void Clear_Selection();

};

template<class T>
class OPENGL_SELECTION_DEBUG_PARTICLES_2D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int index;

    OPENGL_SELECTION_DEBUG_PARTICLES_2D(OPENGL_OBJECT<T>* object)
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::DEBUG_PARTICLES_2D,object) {}

    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

}

#endif
