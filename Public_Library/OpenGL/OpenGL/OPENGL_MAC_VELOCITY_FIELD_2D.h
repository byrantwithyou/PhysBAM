//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_MAC_VELOCITY_FIELD_2D
//#####################################################################
#ifndef __OPENGL_MAC_VELOCITY_FIELD_2D__
#define __OPENGL_MAC_VELOCITY_FIELD_2D__

#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <OpenGL/OpenGL/OPENGL_GRID_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class RANGE;

template<class T_input>
class OPENGL_MAC_VELOCITY_FIELD_2D:public OPENGL_VECTOR_FIELD_2D<ARRAY<VECTOR<T_input,2> > >,public OPENGL_GRID_OBJECT<VECTOR<T_input,2> >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,TV::m> TV_INT;
public:
    using OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >::size;using OPENGL_OBJECT<T>::World_Space_Box;

    enum VELOCITY_MODE { FACE_CENTERED, CELL_CENTERED };
    VELOCITY_MODE velocity_mode;

    GRID<TV> grid;
    ARRAY<T,FACE_INDEX<2> > face_velocities;
    ARRAY_VIEW<T,VECTOR<int,2> > &u,&v;
    ARRAY<TV> vector_field,vector_locations;
    ARRAY<bool,TV_INT> *active_cells;
    ARRAY<bool,FACE_INDEX<TV::m> > *active_faces;
    
    OPENGL_MAC_VELOCITY_FIELD_2D(const GRID<TV> &grid_input);
    virtual ~OPENGL_MAC_VELOCITY_FIELD_2D();

    void Update();  // Call when grid/u/v change
    void Print_Cell_Selection_Info(std::ostream& stream,const TV_INT& cell) const override;
    void Print_Node_Selection_Info(std::ostream& stream,const TV_INT& node) const override;
    void Set_Velocity_Mode(VELOCITY_MODE velocity_mode_input);

    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    // convenience functions
    void Toggle_Velocity_Mode();
};
}
#endif
