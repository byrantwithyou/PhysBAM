//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_COLLECTION
//#####################################################################
#ifndef __FLUID_COLLECTION__
#define __FLUID_COLLECTION__

#include <Grid_Tools/Grids/GRID.h>
#include <Incompressible/Incompressible_Fluids/INCOMPRESSIBLE_FLUID_COLLECTION.h>
namespace PhysBAM{
class VIEWER_DIR;

template<class TV>
class FLUID_COLLECTION
{
    typedef typename TV::SCALAR T;
public:
    const GRID<TV>& grid;
    //COMPRESSIBLE_FLUID_COLLECTION<TV>& compressible_fluid_collection;
    INCOMPRESSIBLE_FLUID_COLLECTION<TV> incompressible_fluid_collection;
    
    FLUID_COLLECTION(const GRID<TV>& grid_input);
    FLUID_COLLECTION(const FLUID_COLLECTION&) = delete;
    void operator=(const FLUID_COLLECTION&) = delete;
    virtual ~FLUID_COLLECTION();

//#####################################################################
    void Write_Output_Files(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const;
    void Read_Output_Files(const VIEWER_DIR& viewer_dir);
    void Initialize_Grids();
//#####################################################################
};
}
#endif
