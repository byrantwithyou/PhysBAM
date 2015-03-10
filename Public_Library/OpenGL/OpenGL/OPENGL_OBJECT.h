//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Geoffrey Irving, Michael Lentine, Neil Molino, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OBJECT
//#####################################################################
#ifndef __OPENGL_OBJECT__
#define __OPENGL_OBJECT__

#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/ROTATION.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <OpenGL/OpenGL/Convert_1d_To_3d.h>
#include <OpenGL/OpenGL/Convert_2d_To_3d.h>
#include <OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h> // just so we get gl in the right order
#include <OpenGL/OpenGL/OPENGL_SLICE.h>
#include <boost/function.hpp>
namespace PhysBAM{

template<class T> class OPENGL_SELECTION;

template<class T>
class OPENGL_OBJECT:public NONCOPYABLE
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
    STREAM_TYPE stream_type;
    HASHTABLE<std::string,OPENGL_CALLBACK> viewer_callbacks;

    OPENGL_OBJECT(STREAM_TYPE stream_type);
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

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer,int buffer_size);
    virtual void Set_Selection(OPENGL_SELECTION<T>* selection);
    virtual void Highlight_Selection(OPENGL_SELECTION<T>* selection);
    virtual void Clear_Highlight();
    virtual void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const;
    virtual RANGE<TV> Selection_Bounding_Box(OPENGL_SELECTION<T>* selection) const;
    virtual OPENGL_SELECTION<T>* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION<T>* old_selection,bool& delete_selection);

    virtual void Set_Slice(OPENGL_SLICE *slice_input);
    virtual void Slice_Has_Changed();

    TV World_Space_Point(const VECTOR<T,1>& object_space_point) const
    {return World_Space_Point(Convert_1d_To_3d(object_space_point));}

    TV World_Space_Point(const VECTOR<T,2>& object_space_point) const
    {return World_Space_Point(Convert_2d_To_3d(object_space_point));}

    TV World_Space_Point(const TV& object_space_point) const
    {return *frame*object_space_point;}

    RANGE<TV> World_Space_Box(const RANGE<VECTOR<T,1> >& object_space_box) const
    {return World_Space_Box(RANGE<TV>(TV(object_space_box.min_corner.x,0,0),TV(object_space_box.max_corner.x,0,0)));}

    RANGE<TV> World_Space_Box(const RANGE<VECTOR<T,2> >& object_space_box) const
    {return World_Space_Box(Convert_2d_To_3d(object_space_box));}

    RANGE<TV> World_Space_Box(const RANGE<TV>& object_space_box) const
    {return ORIENTED_BOX<TV>(object_space_box,*frame).Axis_Aligned_Bounding_Box();}

    void Send_Transform_To_GL_Pipeline() const
    {OpenGL_Translate(frame->t);OpenGL_Rotate(frame->r);}
};
}
#endif
