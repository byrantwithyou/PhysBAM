//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D__
#define __OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T>
class OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;
public:
    GRID<TV> mac_grid;
    ARRAY<TRIPLE<VECTOR<int,2>,VECTOR<T,2>,char> > pseudo_dirichlet_cells;
private:
    T velocity_scale;
    std::string filename;
    int frame_loaded;
    bool valid;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::World_Space_Box;using OPENGL_COMPONENT<T>::stream_type;
    using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D(STREAM_TYPE stream_type,const GRID<TV> &grid,const std::string &filename_input);
    
    bool Valid_Frame(int frame_input) const override;
    bool Is_Up_To_Date(int frame) const override { return valid && frame_loaded == frame; }
    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;
    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();

private:
    void Reinitialize(bool force=false);
};
}
#endif
