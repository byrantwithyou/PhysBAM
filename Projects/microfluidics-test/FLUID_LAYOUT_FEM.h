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

    struct CONNECTION_DATA
    {
        int pid;
        bool collapsed;
    };
    // (vert index,pipe index) -> connection data
    typedef HASHTABLE<PAIR<int,int>,ARRAY<CONNECTION_DATA> > CONNECTION;

    struct ELEMENT_DATA
    {
        int block_id;
    };

    struct BLOCK_DATA
    {
        bool regular;
    };

    struct BC_DATA
    {
        BC_TYPE bc_type;
    };

    TRIANGULATED_AREA<T>& area;
    ARRAY<BLOCK_DATA> blocks;
    ARRAY<ELEMENT_DATA> elem_data;
    ARRAY<BC_DATA> bc;
    HASHTABLE<int,int> bc_map; // particle index -> bc index

    FLUID_LAYOUT_FEM();
    ~FLUID_LAYOUT_FEM();
    void Compute(const PARSE_DATA_FEM<TV>& pd);
    void Generate_End(int i,int pipe,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_Arc(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_Corner(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_Pipe(int pipe,const PARSE_DATA_FEM<TV>& pd,const CONNECTION& con);
    bool Generate_3_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_Triangle_Junction(int i,const VECTOR<int,3>& ends,const VECTOR<int,3>& pipes,
        const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Dump_Mesh() const;
    void Dump_Layout() const;
    void Dump_Input(const PARSE_DATA_FEM<TV>& pd) const;

    ARRAY<int> March_Corner(const TV& start_point,int p1,const ARRAY<int>& side,T unit_length);
    ARRAY<int> March_Arc(int p0,const TV& end_point,const ARRAY<int>& side,const TV& c,T unit_length);
    void Mark_BC(const ARRAY<CONNECTION_DATA>& pindices,BC_TYPE bc_type);
    // return (center, normalized start vec, normalied end vec)
    VECTOR<TV,3> Wedge(const TV& joint,const TV& p0,const TV& p1,int half_width,T unit_length) const;
};

}
#endif
