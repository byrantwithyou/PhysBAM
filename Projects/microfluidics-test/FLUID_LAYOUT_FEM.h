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

    // (vert index,pipe index) -> array of particle id
    typedef HASHTABLE<PAIR<int,int>,ARRAY<int> > CONNECTION;

    struct ELEMENT_DATA
    {
        int block_id;
        int greed;
    };

    struct BLOCK_DATA
    {
        int pipe_id;
    };

    struct BC_DATA
    {
        BC_TYPE bc_type;
        TV bc;
    };

    struct PIPE_DATA
    {
        TV u,v; // u = x1-x0, v=x2-x0
        int first_element;
    };
    
    TRIANGULATED_AREA<T>& area;
    ARRAY<BLOCK_DATA> blocks;
    ARRAY<ELEMENT_DATA> elem_data; // per element
    ARRAY<BC_DATA> bc;
    // unordered (particle index, particle index) -> bc index
    HASHTABLE<PAIR<int,int>,int> bc_map;
    HASHTABLE<int,int> particle_bc_map;
    ARRAY<PIPE_DATA> pipes; // per pipe
    
    ARRAY_VIEW<bool> node_blocks_assigned; // FIXME: I want each new element of node_blocks is initialized with -1 rather than T().
    ARRAY_VIEW<int> node_blocks,pressure_dofs,vel_node_dofs;
    ARRAY<int> edge_blocks,vel_edge_dofs;
    int num_vel_dofs,num_pressure_dofs;

    FLUID_LAYOUT_FEM();
    ~FLUID_LAYOUT_FEM();
    void Compute(const PARSE_DATA_FEM<TV>& pd);
    void Allocate_Dofs();
    void Print_Statistics() const;
    void Generate_End(int i,int pipe,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_2_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con,JOINT_TYPE jt);
    void Generate_3_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_3_Joint_SmallMin(int i,const VECTOR<int,3>& ends,const VECTOR<int,3>& pipes,
        const VECTOR<T,3>& angles,const VECTOR<VECTOR<TV,3>,3>& tri,
        const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_3_Joint_LargeMin(int i,const VECTOR<int,3>& ends,const VECTOR<int,3>& pipes,
        const VECTOR<T,3>& angles,const VECTOR<VECTOR<TV,3>,3>& tri,
        const PARSE_DATA_FEM<TV>& pd,CONNECTION& con);
    void Generate_Pipe(int pipe,const PARSE_DATA_FEM<TV>& pd,const CONNECTION& con);

    void Dump_Mesh() const;
    void Dump_Layout() const;
    void Dump_Input(const PARSE_DATA_FEM<TV>& pd) const;
    void Dump_Edge_Blocks() const;
    void Dump_Node_Blocks() const;
    void Dump_Dofs() const;

    void Mark_BC(const ARRAY<int>& pindices,BC_TYPE bc_type,const TV& value);
    // return (center, normalized start vec, normalied end vec)
    VECTOR<TV,3> Wedge(const TV& joint,const TV& p0,const TV& p1,int half_width,T unit_length) const;
    ARRAY<int> Sample_Interpolated(T s,const ARRAY<int>& side0,const ARRAY<int>& side1,T unit_length);
    void Merge_Interpolated(const ARRAY<int>& left,const ARRAY<int>& right);
    ARRAY<int> Polyline(const ARRAY<TV>& points,T unit_length);
    PAIR<ARRAY<int>,ARRAY<int> > Arc(const TV& c,const TV& p0,const TV& p1,int half_width,T unit_length,
        const TV& dir0,T n0,const TV& dir1,T n1);
    PAIR<ARRAY<int>,ARRAY<int> > Corner(const TV& c,const TV& joint,const TV& p0,const TV& p1,T unit_length,
        const TV& dir0,T n0,const TV& dir1,T n1);
    void Weld(int n,const ARRAY<int>& side0,const ARRAY<int>& side1,T unit_length,ARRAY<int>& f0,ARRAY<int>& f1);
};

}
#endif
