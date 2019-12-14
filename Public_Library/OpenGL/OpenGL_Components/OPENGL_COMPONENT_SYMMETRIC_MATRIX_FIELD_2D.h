//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D__
#define __OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D__

#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_SYMMETRIC_MATRIX_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T>
class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_OBJECT<T>::viewer_callbacks;using OPENGL_COMPONENT<T>::viewer_dir;
    ARRAY<SYMMETRIC_MATRIX<T,2> ,VECTOR<int,2> > field;
    OPENGL_SYMMETRIC_MATRIX_FIELD_2D<T> opengl_symmetric_matrix_field;
    std::string field_filename;
    bool valid;

    OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D(const VIEWER_DIR& viewer_dir,const GRID<TV>& grid,const std::string& field_filename_input)
        :OPENGL_COMPONENT<T>(viewer_dir,"Symmetric Matrix Field 2D"),opengl_symmetric_matrix_field(grid,field),
        field_filename(field_filename_input),valid(false)
    {
        viewer_callbacks.Set("increase_size",{[this](){Increase_Size();},"Increase symmetric matrix size"});
        viewer_callbacks.Set("decrease_size",{[this](){Decrease_Size();},"Decrease symmetric matrix size"});
    }

    virtual ~OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D()
    {}

    void Set_Frame() override
    {Reinitialize();}

    void Set_Draw(bool draw_input=true) override
    {OPENGL_COMPONENT<T>::Set_Draw(draw_input);Reinitialize();}

    void Display() const override
    {if(valid&&draw)opengl_symmetric_matrix_field.Display();}
    
    bool Use_Bounding_Box() const override
    {return draw&&valid;}

    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override
    {if(valid&&draw)return opengl_symmetric_matrix_field.Bounding_Box();return RANGE<VECTOR<T,3> >::Centered_Box();}

    void Increase_Size()
    {opengl_symmetric_matrix_field.size*=(T)1.1;}

    void Decrease_Size()
    {opengl_symmetric_matrix_field.size*=1/(T)1.1;}

//#####################################################################
    void Reinitialize();
//#####################################################################
};
}

#endif
