//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D__
#define __OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D__

#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_SYMMETRIC_MATRIX_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T>
class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;
    ARRAY<SYMMETRIC_MATRIX<T,2> ,VECTOR<int,2> > field;
    OPENGL_SYMMETRIC_MATRIX_FIELD_2D<T> opengl_symmetric_matrix_field;
    std::string field_filename;
    int frame_loaded;
    bool valid;

    OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D(STREAM_TYPE stream_type,const GRID<TV>& grid,const std::string& field_filename_input)
        :OPENGL_COMPONENT<T>(stream_type,"Symmetric Matrix Field 2D"),opengl_symmetric_matrix_field(stream_type,grid,field),
        field_filename(field_filename_input),frame_loaded(-1),valid(false)
    {
        is_animation=FILE_UTILITIES::Is_Animated(field_filename);Reinitialize();
    }

    virtual ~OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D()
    {}

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE
    {return FILE_UTILITIES::Frame_File_Exists(field_filename,frame_input);}

    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE
    {return valid && frame_loaded==frame;}

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE
    {OPENGL_COMPONENT<T>::Set_Frame(frame_input);Reinitialize();}

    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE
    {OPENGL_COMPONENT<T>::Set_Draw(draw_input);Reinitialize();}

    void Display() const PHYSBAM_OVERRIDE
    {if(valid&&draw)opengl_symmetric_matrix_field.Display();}
    
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return draw&&valid;}

    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE
    {if(valid&&draw)return opengl_symmetric_matrix_field.Bounding_Box();return RANGE<VECTOR<T,3> >::Centered_Box();}

    void Increase_Size()
    {opengl_symmetric_matrix_field.size*=(T)1.1;}

    void Decrease_Size()
    {opengl_symmetric_matrix_field.size*=1/(T)1.1;}

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D,Increase_Size,"Increase symmetric matrix size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D,Decrease_Size,"Decrease symmetric matrix size");

//#####################################################################
    void Reinitialize(bool force_load_even_if_not_drawn=false);
//#####################################################################
};
}

#endif
