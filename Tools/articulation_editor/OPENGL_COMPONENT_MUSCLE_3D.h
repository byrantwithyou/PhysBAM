//#####################################################################
// Copyright 2005, Jared Go, Ranjitha Kumar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_MUSCLE_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_MUSCLE_3D__
#define __OPENGL_COMPONENT_MUSCLE_3D__

#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>


namespace PhysBAM
{

template<class T> class ARTICULATION_VISUALIZATION;

template<class T,class RW=T>
class OPENGL_COMPONENT_MUSCLE_3D : public OPENGL_COMPONENT
{
public:
    OPENGL_COMPONENT_MUSCLE_3D(ARTICULATION_VISUALIZATION<T> *owner);
    virtual ~OPENGL_COMPONENT_MUSCLE_3D();
    
    // Owning vis
    ARTICULATION_VISUALIZATION<T> *owner;

    // Should we be hilighed?
    bool selected;
    int targpoint;
    int segment;

    // Calls
    virtual void Display(const int in_color=1) const;
    virtual bool Use_Bounding_Box() const { return false; }
    virtual BOX_3D<float> Bounding_Box() const;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    //virtual void Highlight_Selection(OPENGL_SELECTION *selection) { opengl_triangulated_surface.Highlight_Selection(selection); }
    //virtual void Clear_Highlight() { opengl_triangulated_surface.Clear_Highlight(); }

    void AddPointInFrame(const VECTOR<T,3> &pt, const FRAME_3D<T> &frame, const std::string &bname, const std::string &attached,const std::string &st);    

private:
    
public:
    ARRAY< VECTOR<T,3> > points;
    ARRAY< FRAME_3D<T> > frames;
    ARRAY< std::string > bones;
    ARRAY< std::string > attach;
    ARRAY< std::string > segment_type;
    ARRAY< ARRAY<T> > parameters;
private:    
};

}

#endif
