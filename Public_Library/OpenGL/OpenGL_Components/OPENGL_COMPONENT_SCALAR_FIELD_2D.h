//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SCALAR_FIELD_2D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_SCALAR_FIELD_2D__
#define __OPENGL_COMPONENT_SCALAR_FIELD_2D__

#include <OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>

namespace PhysBAM
{

template<class T,class T2=T>
class OPENGL_COMPONENT_SCALAR_FIELD_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;
    using OPENGL_OBJECT<T>::viewer_callbacks;using OPENGL_COMPONENT<T>::viewer_dir;
    OPENGL_COMPONENT_SCALAR_FIELD_2D(const VIEWER_DIR& viewer_dir,GRID<TV> &grid_input,
        const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,
        const char* info_name,typename OPENGL_SCALAR_FIELD_2D<T,T2>::DRAW_MODE draw_mode_input=OPENGL_SCALAR_FIELD_2D<T,T2>::DRAW_TEXTURE);
    virtual ~OPENGL_COMPONENT_SCALAR_FIELD_2D();

    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    void Print_Selection_Info(std::ostream& output_stream) const override;
    void Toggle_Smooth();
    void Toggle_Draw_Mode();
    void Toggle_Color_Map();

private:
    void Reinitialize();

public:
    OPENGL_SCALAR_FIELD_2D<T,T2>  opengl_scalar_field;

private:
    std::string scalar_field_filename;
    bool valid;
};

}

#endif
