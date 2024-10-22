//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Geoffrey Irving, Michael Lentine, Neil Molino, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OBJECT
//#####################################################################
#ifndef __OPENGL_OBJECT__
#define __OPENGL_OBJECT__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <OpenGL/OpenGL/Convert_1d_To_3d.h>
#include <OpenGL/OpenGL/Convert_2d_To_3d.h>
#include <OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h> // just so we get gl in the right order
#include <OpenGL/OpenGL/OPENGL_SLICE.h>
#include <functional>
namespace PhysBAM{

template<class T>
class OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
public:
    FRAME<TV>* frame; // pointer so you can enslave this body to another's motion
    std::string name;
    bool selectable;
    bool visible;
    bool show_name;
    OPENGL_SLICE *slice;
private:
    FRAME<TV> default_frame;

public:
    HASHTABLE<std::string,OPENGL_CALLBACK> viewer_callbacks;

    OPENGL_OBJECT();
    OPENGL_OBJECT(const OPENGL_OBJECT&) = delete;
    void operator=(const OPENGL_OBJECT&) = delete;
    virtual ~OPENGL_OBJECT();

    void Set_Name(const std::string& name_input)
    {name=name_input;}

    void Enslave_Transform_To(OPENGL_OBJECT& object)
    {frame=object.frame;}

    virtual void Display() const;
    virtual bool Use_Bounding_Box() const;
    virtual RANGE<TV> Bounding_Box() const;
    virtual bool Is_Transparent() const;
    virtual void Turn_Smooth_Shading_On();
    virtual void Turn_Smooth_Shading_Off();

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices);
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers);
    virtual void Clear_Selection();
    virtual void Print_Selection_Info(std::ostream &output_stream) const;
    virtual RANGE<TV> Selection_Bounding_Box() const;
    virtual bool Destroy_Selection_After_Frame_Change();

    virtual void Set_Slice(OPENGL_SLICE *slice_input);
    virtual void Slice_Has_Changed();

    TV World_Space_Point(const TV& object_space_point) const
    {return *frame*object_space_point;}

    TV World_Space_Point(const VECTOR<T,1>& object_space_point) const
    {return *frame*VECTOR<T,3>(object_space_point);}

    TV World_Space_Point(const VECTOR<T,2>& object_space_point) const
    {return *frame*VECTOR<T,3>(object_space_point);}

    RANGE<TV> World_Space_Box(const RANGE<TV>& object_space_box) const
    {return ORIENTED_BOX<TV>(object_space_box,*frame).Axis_Aligned_Bounding_Box();}

    RANGE<TV> World_Space_Box(const RANGE<VECTOR<T,1> >& object_space_box) const
    {return World_Space_Box(Convert_1d_To_3d(object_space_box));}

    RANGE<TV> World_Space_Box(const RANGE<VECTOR<T,2> >& object_space_box) const
    {return World_Space_Box(Convert_2d_To_3d(object_space_box));}

    void Send_Transform_To_GL_Pipeline() const
    {OpenGL_Translate(frame->t);OpenGL_Rotate(frame->r);}
};
}
#endif
