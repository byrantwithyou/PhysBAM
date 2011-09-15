//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef _MAYA_SPLINE_PLUGIN_GEOMETRY
#define _MAYA_SPLINE_PLUGIN_GEOMETRY

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{

class MAYA_SPLINE_PLUGIN_GEOMETRY
{
public:
    PhysBAM::ARRAYS<VECTOR<PhysBAM::VECTOR_3D<float> ,3> > controls;
    PhysBAM::GRID_3D<float> grid;

    MAYA_SPLINE_PLUGIN_GEOMETRY();
    ~MAYA_SPLINE_PLUGIN_GEOMETRY();
    
    void Initialize(const GRID_3D<float>& grid_input,const bool zeroed_input=false);
    int Get_Number_Vertices();
};
//#####################################################################
}
#endif
