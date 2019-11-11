//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_FLUID_COLLECTION
//#####################################################################
#ifndef __COMPRESSIBLE_FLUID_COLLECTION__
#define __COMPRESSIBLE_FLUID_COLLECTION__

#include <Grid_Tools/Grids/GRID.h>
#include <Compressible/Equations_Of_State/EOS.h>
namespace PhysBAM{
class VIEWER_DIR;

template<class TV>
class COMPRESSIBLE_FLUID_COLLECTION
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
public:
    const GRID<TV>& grid;

    EOS<T>* eos;
    ARRAY<bool,TV_INT> psi;
    T_ARRAYS_DIMENSION_SCALAR U;
        
    COMPRESSIBLE_FLUID_COLLECTION(const GRID<TV>& grid_input);
    COMPRESSIBLE_FLUID_COLLECTION(const COMPRESSIBLE_FLUID_COLLECTION&) = delete;
    void operator=(const COMPRESSIBLE_FLUID_COLLECTION&) = delete;
    ~COMPRESSIBLE_FLUID_COLLECTION();

    void Set_Equation_Of_State(EOS<T>* eos_input)
    {eos=eos_input;}
    
    //#####################################################################
    void Write_Output_Files(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const;
    void Read_Output_Files(const VIEWER_DIR& viewer_dir);
    void Initialize_Grids();
    //#####################################################################
};
}
#endif
