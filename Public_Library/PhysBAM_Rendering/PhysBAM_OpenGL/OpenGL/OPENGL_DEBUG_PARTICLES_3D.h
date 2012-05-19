//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_DEBUG_PARTICLES_3D
//#####################################################################
#ifndef __OPENGL_DEBUG_PARTICLES_3D__
#define __OPENGL_DEBUG_PARTICLES_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T>
class OPENGL_DEBUG_PARTICLES_3D:public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
public:
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<DEBUG_OBJECT<TV> >& debug_objects;
    OPENGL_COLOR default_color;
    OPENGL_COLOR velocity_color;
    bool draw_velocities;
    T scale_velocities;
    bool wireframe_only;

    OPENGL_DEBUG_PARTICLES_3D(GEOMETRY_PARTICLES<TV>& particles_input,ARRAY<DEBUG_OBJECT<TV> >& debug_objects_input,const OPENGL_COLOR &color_input = OPENGL_COLOR::White());
    ~OPENGL_DEBUG_PARTICLES_3D();

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;

    OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size) PHYSBAM_OVERRIDE;
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;

    void Select_Point(int index);
    void Select_Points(const ARRAY<int> &indices);
    void Clear_Selection();

};

template<class T>
class OPENGL_SELECTION_DEBUG_PARTICLES_3D:public OPENGL_SELECTION
{
public:
    int index;

    OPENGL_SELECTION_DEBUG_PARTICLES_3D(OPENGL_OBJECT *object) : OPENGL_SELECTION(OPENGL_SELECTION::DEBUG_PARTICLES_3D, object) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
