//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D__
#define __OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_GRID_OBJECT.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T>
class OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D:public OPENGL_COMPONENT<T>,public OPENGL_GRID_OBJECT<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    GRID<TV> mac_grid;
    ARRAY<TRIPLE<VECTOR<int,3>,VECTOR<T,3>,char> > pseudo_dirichlet_cells;
private:
    T velocity_scale;
    std::string filename;
    bool valid;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::viewer_dir;
    using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid,const std::string &filename_input);
    
    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;
    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<TV> Bounding_Box() const override;
    void Print_Cell_Selection_Info(std::ostream& stream,const TV_INT& cell) const override;
    void Print_Node_Selection_Info(std::ostream& stream,const TV_INT& node) const override;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();

private:
    void Reinitialize();
};
}
#endif
