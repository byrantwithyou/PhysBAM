//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_LAYOUT_FEM__
#define __FLUID_LAYOUT_FEM__
#include <Core/Data_Structures/HASHTABLE.h>
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

    struct BC_DATA
    {
        BC_TYPE bc_type;
    };

    std::unique_ptr<TRIANGULATED_AREA<T> > area;
    ARRAY<BLOCK_DATA> blocks; // element index -> block
    HASHTABLE<PAIR<int,int>,BC_DATA> bc; // edge -> bc
    int last_block_id=0;

    FLUID_LAYOUT_FEM();
    ~FLUID_LAYOUT_FEM();
    void Compute(const PARSE_DATA_FEM<TV>& pd);
    void Dump_Mesh() const;
    void Dump_Layout() const;
    void Dump_Input(const PARSE_DATA_FEM<TV>& pd) const;
    // shared_point: (local i, local j) -> particle index; i==-1 means the last edge
    PAIR<ARRAY<int>,ARRAY<int> > Generate_Pipe(const TV& v0,const TV& v1,int half_num_cells,T unit_length,
        const HASHTABLE<PAIR<int,int>,int>& shared_point={});
    ARRAY<int> Weld_Arc(int p0,int p1,const ARRAY<int>& side,const TV& c,T unit_length);
    void Mark_BC(const ARRAY<int>& pindices,BC_TYPE bc_type);
    void Pipe_Joint_Connection(const TV& joint,const TV& p0,const TV& p1,int half_width,T unit_length,
        TV& q0,TV& q1) const;
};

}
#endif
