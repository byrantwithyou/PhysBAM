//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GRID_SETUP__
#define __GRID_SETUP__
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR.h>
#include "COMMON.h"
#include "PARSE_DATA.h"

namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class PARSE_DATA;

template<class TV>
struct FLUID_LAYOUT
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const GRID<TV>& grid;

    struct BLOCK
    {
        int num_dofs;
    };

    struct GRID_DATA
    {
        INDEX_TYPE type;
        int block_id;
        int block_dof;
        T bc_value;
        DOF_ID global_id;
    };

    ARRAY<GRID_DATA,FACE_INDEX<TV::m> > used_faces;
    ARRAY<GRID_DATA,TV_INT> used_cells;
    ARRAY<BLOCK> blocks;
    ARRAY<VECTOR<int,2>,DOF_ID> dof_map;
    int num_vertex_blocks;
    bool quiet;
    
    FLUID_LAYOUT(const GRID<TV>& grid): grid(grid) {}
    
    void Compute(const PARSE_DATA<TV>& pd);
    void Dump_Layout() const;
    void Dump_Dofs() const;
    void Dump_Blocks() const;
    DOF_ID Total_Dofs() const {return dof_map.m;}
    void Allocate_Cross_Section_Blocks_Cells(const RANGE<TV_INT>& box,int dir);
    void Allocate_Cross_Section_Blocks_Faces(const RANGE<TV_INT>& box,int dir);

};

}
#endif
