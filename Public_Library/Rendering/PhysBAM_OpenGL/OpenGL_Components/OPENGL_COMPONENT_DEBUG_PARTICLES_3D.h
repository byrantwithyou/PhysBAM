//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEBUG_PARTICLES_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEBUG_PARTICLES_3D__
#define __OPENGL_COMPONENT_DEBUG_PARTICLES_3D__

#include <Tools/Arrays/ARRAY.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_DEBUG_PARTICLES_3D.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_DEBUG_PARTICLES_3D:public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_COMPONENT_DEBUG_PARTICLES_3D(const std::string &filename);
    virtual ~OPENGL_COMPONENT_DEBUG_PARTICLES_3D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Selection_Bounding_Box(OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    void Show_Colored_Wireframe();
    void Toggle_Draw_Velocities();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Command_Prompt();
    void Set_Slice(OPENGL_SLICE *slice_input);
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_3D,Show_Colored_Wireframe,"Show colored wireframe");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_3D,Toggle_Draw_Velocities,"Toggle draw velocities");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_3D,Increase_Vector_Size,"Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_3D,Decrease_Vector_Size,"Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_3D,Command_Prompt,"Command prompt");

private:
    void Reinitialize(bool force=false);

    void Command_Prompt_Response();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_3D, Command_Prompt_Response, "");

public:
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<DEBUG_OBJECT<TV> >& debug_objects;
    OPENGL_DEBUG_PARTICLES_3D<T>& opengl_particles;

private:
    std::string filename;
    int frame_loaded;
    int set;
    int set_loaded;
    bool valid;
    bool draw_multiple_particle_sets;
};

template<class T>
class OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D : public OPENGL_SELECTION
{
public:
    int index;  // index into particles array
    VECTOR<T,3> location;

    OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D(OPENGL_OBJECT *object) : OPENGL_SELECTION(OPENGL_SELECTION::DEBUG_PARTICLES_3D, object) {}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
