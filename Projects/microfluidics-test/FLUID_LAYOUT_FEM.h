//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_LAYOUT_FEM__
#define __FLUID_LAYOUT_FEM__
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR.h>
#include "PARSE_DATA_FEM.h"

namespace PhysBAM{

template<class T> class TRIANGULATED_AREA;

template<class TV>
struct FLUID_LAYOUT_FEM
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    struct BLOCK_DATA
    {
        int block_id;
        bool regular;
    };

    std::unique_ptr<TRIANGULATED_AREA<T> > area;
    ARRAY<BLOCK_DATA> blocks; // element index -> block

    FLUID_LAYOUT_FEM();
    ~FLUID_LAYOUT_FEM();
    void Compute(const PARSE_DATA_FEM<TV>& pd);
    void Dump_Mesh() const;
    void Dump_Layout() const;
    void Dump_Input(const PARSE_DATA_FEM<TV>& pd) const;
    void Generate_Pipe(const TV& v0,const TV& v1,int half_width,T unit_length);
};

}
#endif
