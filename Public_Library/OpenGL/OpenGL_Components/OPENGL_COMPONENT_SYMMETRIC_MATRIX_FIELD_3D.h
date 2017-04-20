//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D__
#define __OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D__

#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_SYMMETRIC_MATRIX_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T>
class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::is_animation;using OPENGL_COMPONENT<T>::stream_type;
    using OPENGL_OBJECT<T>::viewer_callbacks;
    ARRAY<SYMMETRIC_MATRIX<T,3>,VECTOR<int,3> > field;
    OPENGL_SYMMETRIC_MATRIX_FIELD_3D<T> opengl_symmetric_matrix_field;
    std::string field_filename;
    int frame_loaded;
    bool valid;

    OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D(STREAM_TYPE stream_type,const GRID<TV>& grid,const std::string& field_filename_input)
        :OPENGL_COMPONENT<T>(stream_type,"Symmetric Matrix Field 3D"),opengl_symmetric_matrix_field(stream_type,grid,field),
        field_filename(field_filename_input),frame_loaded(-1),valid(false)
    {
        is_animation=Is_Animated(field_filename);Reinitialize();

        viewer_callbacks.Set("increase_size",{[this](){Increase_Size();},"Increase symmetric matrix size"});
        viewer_callbacks.Set("decrease_size",{[this](){Decrease_Size();},"Decrease symmetric matrix size"});
    }

    virtual ~OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D()
    {}

    bool Valid_Frame(int frame_input) const override
    {return Frame_File_Exists(field_filename,frame_input);}

    bool Is_Up_To_Date(int frame) const override
    {return valid && frame_loaded==frame;}

    void Set_Frame(int frame_input) override
    {OPENGL_COMPONENT<T>::Set_Frame(frame_input);Reinitialize();}

    void Set_Draw(bool draw_input=true) override
    {OPENGL_COMPONENT<T>::Set_Draw(draw_input);Reinitialize();}

    void Display() const override
    {if(valid&&draw)opengl_symmetric_matrix_field.Display();}
    
    bool Use_Bounding_Box() const override
    {return draw&&valid;}

    virtual RANGE<TV> Bounding_Box() const override
    {if(valid&&draw)return opengl_symmetric_matrix_field.Bounding_Box();return RANGE<VECTOR<T,3> >::Centered_Box();}

    void Set_Slice(OPENGL_SLICE* slice_input) override
    {slice=slice_input;opengl_symmetric_matrix_field.Set_Slice(slice_input);}

    void Slice_Has_Changed() override
    {opengl_symmetric_matrix_field.Slice_Has_Changed();} 

    void Increase_Size()
    {opengl_symmetric_matrix_field.size*=(T)1.1;}

    void Decrease_Size()
    {opengl_symmetric_matrix_field.size*=1/(T)1.1;}

//#####################################################################
    void Reinitialize(bool force_load_even_if_not_drawn=false);
//#####################################################################
};
}

#endif
