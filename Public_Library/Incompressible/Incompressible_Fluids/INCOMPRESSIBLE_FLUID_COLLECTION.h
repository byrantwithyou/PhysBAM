//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_FLUID_COLLECTION
//#####################################################################
#ifndef __INCOMPRESSIBLE_FLUID_COLLECTION__
#define __INCOMPRESSIBLE_FLUID_COLLECTION__

#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Incompressible/Grid_Based_Fields/DENSITY_CONTAINER.h>
#include <Incompressible/Grid_Based_Fields/TEMPERATURE_CONTAINER.h>
namespace PhysBAM{
class VIEWER_DIR;

template<class TV>
class INCOMPRESSIBLE_FLUID_COLLECTION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    const GRID<TV>& grid;
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities;
    ARRAY<T,TV_INT> viscosity;
    //DENSITY_CONTAINER<TV> density_container;
    //TEMPERATURE_CONTAINER<TV> temperature_container;
    
    INCOMPRESSIBLE_FLUID_COLLECTION(const GRID<TV>& grid_input);
    INCOMPRESSIBLE_FLUID_COLLECTION(const INCOMPRESSIBLE_FLUID_COLLECTION&) = delete;
    void operator=(const INCOMPRESSIBLE_FLUID_COLLECTION&) = delete;
    ~INCOMPRESSIBLE_FLUID_COLLECTION();

//#####################################################################
    void Write_Output_Files(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const;
    void Read_Output_Files(const VIEWER_DIR& viewer_dir);
    void Initialize_Grids();
//#####################################################################
};
}
#endif
