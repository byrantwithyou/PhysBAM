//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D
//#####################################################################
#ifndef __OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D__
#define __OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
namespace PhysBAM{

template<class T>
class OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    using OPENGL_OBJECT<T>::World_Space_Box;
    T height_scale;
    GRID<TV>& grid;
    ARRAY<TV,TV_INT >& V_minus;
    ARRAY<TV,TV_INT >& V_plus;
    LEVELSET<TV>& levelset;
    OPENGL_VECTOR_FIELD_3D<T> minus;
    OPENGL_VECTOR_FIELD_3D<T> plus;

    OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D(STREAM_TYPE stream_type,GRID<TV>& grid,ARRAY<TV,TV_INT>& V_minus,ARRAY<TV,TV_INT>& V_plus,LEVELSET<TV>& levelset);
    virtual ~OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D();

//#####################################################################
    void Update();  // Call when grid/V change
    void Display() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    // convenience functions
    void Scale_Vector_Size(const T scale);
    void Scale_Height(const T height_scale);
//#####################################################################
};
}
#endif
