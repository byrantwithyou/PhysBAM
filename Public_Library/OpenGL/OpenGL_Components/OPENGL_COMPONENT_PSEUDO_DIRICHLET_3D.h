//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D__
#define __OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T>
class OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;
public:
    GRID<TV> mac_grid;
    ARRAY<TRIPLE<VECTOR<int,3>,VECTOR<T,3>,char> > pseudo_dirichlet_cells;
private:
    T velocity_scale;
    std::string filename;
    int frame_loaded;
    bool valid;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D(STREAM_TYPE stream_type,const GRID<TV> &grid,const std::string &filename_input);
    
    bool Valid_Frame(int frame_input) const override;
    bool Is_Up_To_Date(int frame) const override { return valid && frame_loaded == frame; }
    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;
    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* current_selection) const override;

    void Set_Vector_Size(const T vector_size);

    void Increase_Vector_Size();
    void Decrease_Vector_Size();

private:
    void Reinitialize(bool force=false);
};
}
#endif
