//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_FLUID_COLLECTION
//#####################################################################
#ifndef __INCOMPRESSIBLE_FLUID_COLLECTION__
#define __INCOMPRESSIBLE_FLUID_COLLECTION__

#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Parallel_Computation/THREADED_UNIFORM_GRID.h>
#include <Incompressible/Grid_Based_Fields/DENSITY_CONTAINER.h>
#include <Incompressible/Grid_Based_Fields/TEMPERATURE_CONTAINER.h>
namespace PhysBAM{

template<class TV>
class INCOMPRESSIBLE_FLUID_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
public:
    const GRID<TV>& grid;
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities;
    ARRAY<T,TV_INT> viscosity;
    //DENSITY_CONTAINER<TV> density_container;
    //TEMPERATURE_CONTAINER<TV> temperature_container;
    
    INCOMPRESSIBLE_FLUID_COLLECTION(const GRID<TV>& grid_input);
    ~INCOMPRESSIBLE_FLUID_COLLECTION();

//#####################################################################
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const;
    void Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame);
    void Initialize_Grids();
    void Sync_Data(INCOMPRESSIBLE_FLUID_COLLECTION<TV>& fluid_collection,THREADED_UNIFORM_GRID<TV>& threaded_grid);
    void Distribute_Data(INCOMPRESSIBLE_FLUID_COLLECTION<TV>& fluid_collection,THREADED_UNIFORM_GRID<TV>& threaded_grid);
//#####################################################################
};
}
#endif
