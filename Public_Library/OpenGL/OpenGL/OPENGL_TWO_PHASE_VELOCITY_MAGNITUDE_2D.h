//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D
//#####################################################################
#ifndef __OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D__
#define __OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
namespace PhysBAM{

template<class T>
class OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D:public OPENGL_OBJECT
{
    typedef VECTOR<T,2> TV;
public:
    T height_scale;
    GRID<TV>& grid;
    ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_minus;
    ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_plus;
    LEVELSET<TV>& levelset;
    OPENGL_VECTOR_FIELD_3D<T> minus;
    OPENGL_VECTOR_FIELD_3D<T> plus;

    OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D(GRID<TV>& grid,ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_minus,ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_plus,LEVELSET<TV>& levelset);
    virtual ~OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D();

//#####################################################################
    void Update();  // Call when grid/V change
    void Display() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    // convenience functions
    void Scale_Vector_Size(const T scale);
    void Scale_Height(const T height_scale);
//#####################################################################
};
}
#endif
