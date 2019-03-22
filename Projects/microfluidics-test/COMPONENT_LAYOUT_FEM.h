//#####################################################################
// Copyright 2018.
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
#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_VECTOR.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <map>
#include "BLOCK_MATRIX.h"
#include "BLOCK_VECTOR.h"
#include "CANONICAL_COMPONENT.h"
#include "COMMON.h"
#include "COMPONENT_PIPE.h"
#include "XFORM.h"

namespace PhysBAM{

template<class TV> struct COMPONENT_LAYOUT_FEM;
template<class T> struct CACHED_ELIMINATION_MATRIX;
template<class TV> struct BLOCK_MESHING_ITERATOR;

PHYSBAM_DECLARE_ELEMENT_ID(COMPONENT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(SEPARATOR_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REFERENCE_BLOCK_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REFERENCE_CONNECTION_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REFERENCE_IRREGULAR_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(IRREG_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);

extern double comp_tol;

template<class T>
struct COMPONENT_LAYOUT_FEM<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<int,3> IV3;
    typedef VECTOR<int,2> IV;

    T target_length;
    T mu=1;
    T unit_m,unit_s,unit_kg;
    
    // BLOCKS

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

    struct BLOCK
    {
        CANONICAL_BLOCK<T>* block;
        XFORM<TV> xform;
        ARRAY<BLOCK_CONNECTION,CON_ID> connections;
        ARRAY<IRREG_ID> edge_on; // for edge-on (index in irregular_connections)
        int flags=0; // 1=separator, 2=separator-eligible
        REFERENCE_BLOCK_ID ref_id=REFERENCE_BLOCK_ID(-7);
        DOF_COUNTS num_dofs;
    };
    ARRAY<BLOCK,BLOCK_ID> blocks;
    CANONICAL_BLOCK<T>* empty_canonical_block;

    // IRREGULAR CONNECTIONS

    // regular is master
    struct IRREGULAR_CONNECTION
    {
        BLOCK_ID regular=BLOCK_ID(-7);
        CON_ID con_id=CON_ID(-7);
        // one for each dof on cross section, starting from owned side of cross section
        ARRAY<PAIR<BLOCK_ID,int> > edge_on_v,edge_on_e;
        REFERENCE_IRREGULAR_ID ref_id=REFERENCE_IRREGULAR_ID(-7);
    };

    ARRAY<IRREGULAR_CONNECTION,IRREG_ID> irregular_connections;

    // BOUNDARY CONDITIONS

    struct BOUNDARY_CONDITION
    {
        BLOCK_ID b=BLOCK_ID(-7);
        INTERVAL<int> bc_v,bc_e;
        ARRAY<TV> data_v,data_e;
        TV normal;
    };
    ARRAY<BOUNDARY_CONDITION> bc_v,bc_t;

    ANALYTIC_VECTOR<TV>* analytic_velocity=0;
    ANALYTIC_SCALAR<TV>* analytic_pressure=0;

    TV Traction(const TV& N,const TV& X) const
    {
        SYMMETRIC_MATRIX<T,TV::m> stress=analytic_velocity->dX(X/unit_m,0).Twice_Symmetric_Part()/unit_s*mu;
        stress-=analytic_pressure->f(X/unit_m,0)*unit_kg/sqr(unit_s);
        return stress*N;
    }

    TV Force(const TV& X) const
    {
        SYMMETRIC_TENSOR<T,0,TV::m> ddU=analytic_velocity->ddX(X/unit_m,0)/(unit_m*unit_s);
        TV f=analytic_pressure->dX(X/unit_m,0)*unit_kg/(unit_m*sqr(unit_s));
        f-=mu*(Contract<1,2>(ddU)+Contract<0,2>(ddU));
        return f;
    }
    
    // DOF MAPPING
    struct DOF_PAIRS
    {
        DOF_COUNTS num_dofs_d,num_dofs_s;
        ARRAY<IV> v,e,p;
    };

    struct REFERENCE_BLOCK_DATA
    {
        BLOCK_ID b;
        DOF_COUNTS num_dofs;
        ARRAY<int> dof_map_v,dof_map_e,dof_map_p;
        DOF_PAIRS pairs;
        ARRAY<DOF_PAIRS,CON_ID> regular_pairs; // this block is the source, the connection is the destination
        int mat_id=-7;
    };

    ARRAY<REFERENCE_BLOCK_DATA,REFERENCE_BLOCK_ID> reference_block_data;

    typedef std::tuple<REFERENCE_BLOCK_ID,CON_ID,REFERENCE_BLOCK_ID,CON_ID> REGULAR_CON_KEY;
    HASHTABLE<REGULAR_CON_KEY,REFERENCE_CONNECTION_ID> regular_connection_hash;

    struct REFERENCE_CONNECTION_DATA
    {
        BLOCK_ID b[2];
        CON_ID con_id[2];
        int mat_id=-7;
    };

    ARRAY<REFERENCE_CONNECTION_DATA,REFERENCE_CONNECTION_ID> reference_connection_data;

    struct REFERENCE_IRREGULAR_DATA_HELPER
    {
        BLOCK_ID b=BLOCK_ID(-7);
        DOF_PAIRS irreg_pairs[2]; // dof[to]; 0=regular, 1=irregular; from=1-to
        int mat_id=-7;
    };

    struct REFERENCE_IRREGULAR_DATA
    {
        IRREG_ID ic_id;
        ARRAY<REFERENCE_IRREGULAR_DATA_HELPER> pairs;
    };

    ARRAY<REFERENCE_IRREGULAR_DATA,REFERENCE_IRREGULAR_ID> reference_irregular_data;

    void Compute();
    void Update_Masters();

    // return: flags indicating which connections interact
    // pair: block + master mask
    HASHTABLE<PAIR<CANONICAL_BLOCK<T>*,int>,int> separates_dofs;
    int Separates_Dofs(BLOCK_ID b);
    HASHTABLE<std::tuple<CANONICAL_BLOCK<T>*,CON_ID,CANONICAL_BLOCK<T>*,CON_ID>,PAIR<CANONICAL_BLOCK<T>*,ARRAY<int> > > merge_canonical_blocks;
    void Merge_Blocks(BLOCK_ID id,CON_ID con_id);
    PAIR<CANONICAL_BLOCK<T>*,ARRAY<int> >*
        Merge_Canonical_Blocks(CANONICAL_BLOCK<T>* id0,CON_ID con_id0,
            XFORM<TV> xf0,CANONICAL_BLOCK<T>* id1,CON_ID con_id1,XFORM<TV> xf1);
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

    // MATRIX ASSEMBLY

    HASHTABLE<CANONICAL_BLOCK<T>*,BLOCK_MATRIX<T> > canonical_block_matrices;
    ARRAY<BLOCK_VECTOR<T>,BLOCK_ID> rhs_block_list;
    ARRAY<TRIPLE<BLOCK_ID,BLOCK_ID,int> > nonzero_blocks;

    void Compute_Matrix_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Fill_Canonical_Block_Matrix(BLOCK_MATRIX<T>& mat,const CANONICAL_BLOCK<T>* cb);
    void Fill_Block_Matrix(BLOCK_MATRIX<T>& M,const REFERENCE_BLOCK_DATA& rd);
    void Fill_Connection_Matrix(BLOCK_MATRIX<T>& M,const REFERENCE_CONNECTION_DATA& cd);
    void Fill_Irregular_Connection_Matrix(ARRAY<BLOCK_MATRIX<T> >& M,const REFERENCE_IRREGULAR_DATA& ri);
    void Copy_Matrix_Data(BLOCK_MATRIX<T>& A,BLOCK_ID b,
        const DOF_PAIRS& dpa,const DOF_PAIRS& dpb,BLOCK_ID ar,BLOCK_ID ac) const;
    void Copy_Vector_Data(const BLOCK_VECTOR<T>& B,BLOCK_ID b,const DOF_PAIRS& dp);
    void Init_Block_Matrix(BLOCK_MATRIX<T>& M,BLOCK_ID a,BLOCK_ID b) const;
    void Init_Block_Vector(BLOCK_VECTOR<T>& M,BLOCK_ID b) const;
    void Init_Block_Vector(BLOCK_VECTOR<T>& M,const CANONICAL_BLOCK<T>* cb) const;
    void Times_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& v,const BLOCK_VECTOR<T>& u) const;
    void Times_P_U(BLOCK_ID b,BLOCK_VECTOR<T>& w,const ARRAY<T>& div_v,const ARRAY<T>& div_e) const;
    void Times_Line_Integral_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& w,const BLOCK_VECTOR<T>& u) const;
    void Apply_To_RHS(BLOCK_ID b,const BLOCK_VECTOR<T>& w);


    // ELIMINATION

    void Eliminate_Irregular_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Eliminate_Non_Seperators(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Eliminate_Strip(CACHED_ELIMINATION_MATRIX<T>& cem,const ARRAY<BLOCK_ID>& a);
    void Eliminate_Simple(CACHED_ELIMINATION_MATRIX<T>& cem,BLOCK_ID first,CON_ID con_id_source);


    // PARSING AND ASSEMBLY

    void Parse_Input(const std::string& pipe_file);
    MATRIX<T,2> Compute_Xform(const TV& dir); // dir is normalized

    struct VERTEX_DATA
    {
        TV X;
        BLOCK_CONNECTION con;
    };

    void Emit_Component_Blocks(const CANONICAL_COMPONENT<T>* cc,const XFORM<TV>& xf,ARRAY<VERTEX_DATA>& vd);
    void Set_Connector(VERTEX_DATA& vd,BLOCK_ID id,CON_ID con_id);


    // DEBUGGING

    void Visualize_Block_State(BLOCK_ID b) const;
    void Transform_Solution(const CACHED_ELIMINATION_MATRIX<T>& cem,bool transpose);
    void Visualize_Solution(BLOCK_ID b) const;
    void Dump_World_Space_System(const CACHED_ELIMINATION_MATRIX<T>& cem) const;
    void Dump_World_Space_Vector(const char* name) const;
    void Transform_To_World_Space(BLOCK_MATRIX<T>& M,const BLOCK_MATRIX<T>& B,BLOCK_ID a,BLOCK_ID b) const;
    void Visualize_Flat_Dofs() const;

    // OTHER

    ~COMPONENT_LAYOUT_FEM();
    RANGE<TV> Compute_Bounding_Box() const;
};

}
#endif
