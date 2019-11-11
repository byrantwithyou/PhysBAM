//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPONENT_LAYOUT_FEM__
#define __COMPONENT_LAYOUT_FEM__
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include "CANONICAL_COMPONENT.h"
#include "COMMON.h"
#include "COMPONENT_PIPE.h"
#include "XFORM.h"
#include <map>

namespace PhysBAM{

template<class TV> struct BLOCK_MESHING_ITERATOR;

PHYSBAM_DECLARE_ELEMENT_ID(COMPONENT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(SEPARATOR_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REFERENCE_BLOCK_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REFERENCE_CONNECTION_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REFERENCE_IRREGULAR_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(IRREG_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(RID_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);

extern double comp_tol;

struct BLOCK_CONNECTION
{
    BLOCK_ID id;
    CON_ID con_id;
    IRREG_ID irreg_id;
    bool master,is_regular;

    BLOCK_CONNECTION(BLOCK_ID b=BLOCK_ID(-7),CON_ID c=CON_ID(-7))
        :id(b),con_id(c),irreg_id(-7),master(false),is_regular(true)
    {}

    BLOCK_CONNECTION(IRREG_ID i)
        :id(-7),con_id(-7),irreg_id(i),master(false),is_regular(false)
    {}

    void Set_Irreg(IRREG_ID i)
    {
        con_id=CON_ID(-7);
        irreg_id=i;
        is_regular=false;
    }
};

template<class T>
struct BLOCK
{
    typedef VECTOR<T,2> TV;
    CANONICAL_BLOCK<T>* block;
    XFORM<TV> xform;
    ARRAY<BLOCK_CONNECTION,CON_ID> connections;
    ARRAY<PAIR<IRREG_ID,int> > edge_on; // for edge-on (index in irregular_connections and edge_on)
    int flags=0; // 1=separator, 2=separator-eligible
    REFERENCE_BLOCK_ID ref_id=REFERENCE_BLOCK_ID(-7);
};

struct IRREGULAR_EDGE_DATA
{
    BLOCK_ID b;
    int e,v0,v1; // dofs; v0 borders with previous array entry
};

// regular is master
struct IRREGULAR_CONNECTION
{
    BLOCK_ID regular=BLOCK_ID(-7);
    CON_ID con_id=CON_ID(-7);
    // one for each edge on cross section, starting from side owned by regular block
    ARRAY<IRREGULAR_EDGE_DATA> edge_on;
    REFERENCE_IRREGULAR_ID ref_id=REFERENCE_IRREGULAR_ID(-7);
};

template<class T>
struct BOUNDARY_CONDITION
{
    typedef VECTOR<T,2> TV;
    BLOCK_ID b=BLOCK_ID(-7);
    INTERVAL<int> bc_v,bc_e;
    T flow_rate;
    TV traction;
    TV normal;
};

struct DOF_PAIRS
{
    typedef VECTOR<int,2> IV;
    DOF_COUNTS num_dofs_d,num_dofs_s;
    ARRAY<IV> v,e,p;
};

struct REFERENCE_BLOCK_DATA
{
    BLOCK_ID b;
    DOF_COUNTS num_dofs_d,num_dofs_s;
    ARRAY<int> dof_map_v,dof_map_e,dof_map_p;
    DOF_PAIRS pairs;
    ARRAY<int> ticks_e; // edge index e -> k. tick is on the side of S(e)(k)
    ARRAY<int> ticks_t; // triangle tick masks
};

struct REFERENCE_CONNECTION_DATA
{
    BLOCK_ID b[2];
    CON_ID con_id[2];
    DOF_PAIRS reg_pairs[2]; // dof[dest-index], from-index=1-(dest-index)
};

struct REFERENCE_IRREGULAR_DATA_HELPER
{
    BLOCK_ID b=BLOCK_ID(-7);
    DOF_PAIRS irreg_pairs[2]; // dof[to]; 0=regular, 1=irregular; from=1-to
};

struct REFERENCE_IRREGULAR_DATA
{
    IRREG_ID ic_id;
    ARRAY<REFERENCE_IRREGULAR_DATA_HELPER,RID_ID> pairs;
    ARRAY<PAIR<RID_ID,bool> > mapping;
};

template<class T>
struct COMPONENT_LAYOUT_FEM
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<int,3> IV3;
    typedef VECTOR<int,2> IV;

    T target_length;
    T depth;
    int depth_layers;
    T unit_m,unit_s,unit_kg;
    bool force_blk_ref=false;

    int num_merge=0;
    
    // BLOCKS
    ARRAY<BLOCK<T>,BLOCK_ID> blocks;
    CANONICAL_BLOCK<T>* empty_canonical_block;

    // IRREGULAR CONNECTIONS

    ARRAY<IRREGULAR_CONNECTION,IRREG_ID> irregular_connections;

    // BOUNDARY CONDITIONS

    ARRAY<BOUNDARY_CONDITION<T> > bc_v,bc_t;
    
    // DOF MAPPING

    ARRAY<REFERENCE_BLOCK_DATA,REFERENCE_BLOCK_ID> reference_block_data;

    typedef std::tuple<REFERENCE_BLOCK_ID,CON_ID,REFERENCE_BLOCK_ID,CON_ID,ARRAY<REFERENCE_IRREGULAR_ID> > REGULAR_CON_KEY;
    HASHTABLE<REGULAR_CON_KEY,REFERENCE_CONNECTION_ID> regular_connection_hash;

    ARRAY<REFERENCE_CONNECTION_DATA,REFERENCE_CONNECTION_ID> reference_connection_data;

    ARRAY<REFERENCE_IRREGULAR_DATA,REFERENCE_IRREGULAR_ID> reference_irregular_data;

    void Update_Masters(); // L

    // return: flags indicating which connections interact
    // pair: block + master mask
    // L
    HASHTABLE<PAIR<CANONICAL_BLOCK<T>*,int>,int> separates_dofs;
    int Separates_Dofs(BLOCK_ID b);
    HASHTABLE<std::tuple<CANONICAL_BLOCK<T>*,CON_ID,CANONICAL_BLOCK<T>*,CON_ID>,TRIPLE<CANONICAL_BLOCK<T>*,ARRAY<int>,ARRAY<int> > > merge_canonical_blocks;
    template<class F>
    void Merge_Blocks(BLOCK_ID id,CON_ID con_id,BLOCK_ID id2,F alias);
    TRIPLE<CANONICAL_BLOCK<T>*,ARRAY<int>,ARRAY<int> >& Merge_Canonical_Blocks(
        CANONICAL_BLOCK<T>* cb0,CON_ID con_id0,XFORM<TV> xf0,
        CANONICAL_BLOCK<T>* cb1,CON_ID con_id1,XFORM<TV> xf1);
    int Approx_Dof_Count(BLOCK_ID b);
    void Merge_Blocks();
    void Compute_Reference_Blocks();
    void Compute_Reference_Regular_Connections();
    void Compute_Reference_Irregular_Connections();
    int Compute_Connection_Hash(BLOCK_ID b0,CON_ID con_id0,BLOCK_ID b1,CON_ID con_id1);
    PAIR<int,int> Remap_Owned_Dofs(ARRAY<int>& map_v,ARRAY<int>& map_e,BLOCK_ID b);
    void Compute_Dof_Remapping(REFERENCE_BLOCK_DATA& rd);
    void Compute_Dof_Pairs();
    void Compute_Dof_Pairs(REFERENCE_BLOCK_DATA& rd);
    void Compute_Dof_Pairs(REFERENCE_CONNECTION_DATA& rc);
    void Compute_Dof_Pairs(REFERENCE_IRREGULAR_DATA& ri);
    void Fill_Num_Dofs(DOF_PAIRS& dp,BLOCK_ID d,BLOCK_ID s);
    // b0 side is master
    REGULAR_CON_KEY Regular_Connection_Key(BLOCK_ID b0,CON_ID c0,BLOCK_ID b1) const;
    const DOF_PAIRS& Regular_Connection_Pair(BLOCK_ID b,CON_ID con_id,bool is_dest) const;
    void Fill_Element_Tick_Masks();
    void Fill_Reference_Ticks();

    // OTHER

    ~COMPONENT_LAYOUT_FEM();
    RANGE<TV> Compute_Bounding_Box() const;
};

}
#endif
