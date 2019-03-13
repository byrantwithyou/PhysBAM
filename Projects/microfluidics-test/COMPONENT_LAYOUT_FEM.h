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
#include "COMMON.h"
#include "XFORM.h"

namespace PhysBAM{

template<class TV> struct COMPONENT_LAYOUT_FEM;
template<class T> struct CACHED_ELIMINATION_MATRIX;
template<class TV> struct BLOCK_MESHING_ITERATOR;

PHYSBAM_DECLARE_ELEMENT_ID(COMPONENT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(SEPARATOR_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(CANONICAL_BLOCK_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REFERENCE_BLOCK_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REFERENCE_CONNECTION_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REFERENCE_IRREGULAR_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(IRREG_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(CON_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);

template<class T>
struct COMPONENT_LAYOUT_FEM<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<int,3> IV3;
    typedef VECTOR<int,2> IV;

    static constexpr T comp_tol=(T)1e-10;
    T target_length;
    T mu=1;
    
    // v = range of vertex indices in cross section
    // e = range of edge indices in cross section
    // if own_first, first half of v and e is owned by this block
    // if v or e has odd size; middle is master
    struct CROSS_SECTION
    {
        INTERVAL<int> v,e;
        bool own_first;
    };

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
        SYMMETRIC_MATRIX<T,TV::m> stress=analytic_velocity->dX(X,0).Twice_Symmetric_Part()*mu;
        stress-=analytic_pressure->f(X,0);
        return stress*N;
    }

    TV Force(const TV& X) const
    {
        SYMMETRIC_TENSOR<T,0,TV::m> ddU=analytic_velocity->ddX(X,0);
        TV f=analytic_pressure->dX(X,0);
        f-=mu*(Contract<1,2>(ddU)+Contract<0,2>(ddU));
        return f;
    }

    struct CANONICAL_BLOCK
    {
        ARRAY<CROSS_SECTION,CON_ID> cross_sections;
        ARRAY<TV> X;
        ARRAY<IV3> E;
        ARRAY<IV> S;
        ARRAY<int> bc_v,bc_e;
    };
    ARRAY<CANONICAL_BLOCK,CANONICAL_BLOCK_ID> canonical_blocks;

    ARRAY<BLOCK_MATRIX<T>,CANONICAL_BLOCK_ID> canonical_block_matrices;

    ARRAY<BLOCK_VECTOR<T>,BLOCK_ID> rhs_block_list;

    ARRAY<TRIPLE<BLOCK_ID,BLOCK_ID,int> > nonzero_blocks;
    
    // first is master, second is slave
    struct BLOCK_CONNECTION
    {
        BLOCK_ID id;
        CON_ID con_id;
        IRREG_ID irreg_id;
        bool master,is_regular;

        BLOCK_CONNECTION(BLOCK_ID b=BLOCK_ID(-7),CON_ID c=CON_ID(-7))
            :id(b),con_id(c),irreg_id(-7),master(false),is_regular(true)
        {}

        BLOCK_CONNECTION(BLOCK_ID b,IRREG_ID i)
            :id(b),con_id(-7),irreg_id(i),master(false),is_regular(false)
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
        CANONICAL_BLOCK_ID block;
        XFORM<TV> xform;
        ARRAY<BLOCK_CONNECTION,CON_ID> connections;
        ARRAY<IRREG_ID> edge_on; // for edge-on (index in irregular_connections)
        int flags=0; // 1=separator, 2=separator-eligible
        REFERENCE_BLOCK_ID ref_id=REFERENCE_BLOCK_ID(-7);
    };

    // regular is master
    struct IRREGULAR_CONNECTION
    {
        BLOCK_ID regular=BLOCK_ID(-7);
        CON_ID con_id;
        // one for each dof on cross section, starting from owned side of cross section
        ARRAY<PAIR<BLOCK_ID,int> > edge_on_v,edge_on_e;
        IRREG_ID ref_ic=IRREG_ID(-7);
        REFERENCE_IRREGULAR_ID ref_id=REFERENCE_IRREGULAR_ID(-7);
    };

    // neighbor block i is given index ~i and con_id=-1.
    struct CANONICAL_COMPONENT
    {
        ARRAY<BLOCK,BLOCK_ID> blocks;
        ARRAY<IRREGULAR_CONNECTION,IRREG_ID> irregular_connections;
    };

    ARRAY<BLOCK,BLOCK_ID> blocks;
    ARRAY<IRREGULAR_CONNECTION,IRREG_ID> irregular_connections;

    struct PIPE_KEY
    {
        int num_dofs;
        T width;
        T length;
        bool operator<(const PIPE_KEY& p) const
        {
            if(num_dofs<p.num_dofs) return true;
            if(p.num_dofs<num_dofs) return false;
            if(width<p.width-comp_tol) return true;
            if(p.width<width-comp_tol) return false;
            if(length<p.length-comp_tol) return true;
            return false;
        }
    };
    std::map<PIPE_KEY,CANONICAL_COMPONENT*> canonical_pipes;
    std::map<PIPE_KEY,CANONICAL_BLOCK_ID> canonical_pipe_blocks;

    struct JOINT_KEY
    {
        int num_dofs;
        T width;
        ARRAY<T> angles;
        bool operator<(const JOINT_KEY& p) const
        {
            if(num_dofs<p.num_dofs) return true;
            if(p.num_dofs<num_dofs) return false;
            if(width<p.width-comp_tol) return true;
            if(p.width<width-comp_tol) return false;
            if(angles.m!=p.angles.m) return angles.m<p.angles.m;
            for(int i=0;i<angles.m;i++)
            {
                if(angles(i)<p.angles(i)-comp_tol) return true;
                if(p.angles(i)<angles(i)-comp_tol) return false;
            }
            return false;
        }
    };
    std::map<JOINT_KEY,PAIR<CANONICAL_COMPONENT*,ARRAY<T> > > canonical_joints;

    struct PIPE_CHANGE_KEY
    {
        int num_dofs[2];
        T width[2];
        T length;

        bool operator<(const PIPE_CHANGE_KEY& p) const
        {
            for(int i=0;i<2;i++)
            {
                if(num_dofs[i]<p.num_dofs[i]) return true;
                if(p.num_dofs[i]<num_dofs[i]) return false;
                if(width[i]<p.width[i]-comp_tol) return true;
                if(p.width[i]<width[i]-comp_tol) return false;
            }
            if(length<p.length-comp_tol) return true;
            return false;
        }
    };
    std::map<PIPE_CHANGE_KEY,CANONICAL_COMPONENT*> canonical_changes;
    std::map<PIPE_CHANGE_KEY,CANONICAL_BLOCK_ID> canonical_change_blocks;

    struct DOF_PAIRS
    {
        ARRAY<IV> v,e,p;
    };

    struct REFERENCE_BLOCK_DATA
    {
        BLOCK_ID b;
        int num_dofs_v=-7,num_dofs_e=-7,num_dofs_p=-7;
        ARRAY<int> dof_map_v,dof_map_e,dof_map_p;
        DOF_PAIRS pairs;
        ARRAY<DOF_PAIRS,CON_ID> regular_pairs; // this block is the source, the connection is the destination
        BLOCK_MATRIX<T> M;
        int mat_id=-7;
    };

    ARRAY<REFERENCE_BLOCK_DATA,REFERENCE_BLOCK_ID> reference_block_data;

    typedef std::tuple<REFERENCE_BLOCK_ID,CON_ID,REFERENCE_BLOCK_ID,CON_ID> REGULAR_CON_KEY;
    HASHTABLE<REGULAR_CON_KEY,REFERENCE_CONNECTION_ID> regular_connection_hash;

    struct REFERENCE_CONNECTION_DATA
    {
        BLOCK_ID b[2];
        CON_ID con_id[2];

        BLOCK_MATRIX<T> M;
        int mat_id=-7;
    };

    ARRAY<REFERENCE_CONNECTION_DATA,REFERENCE_CONNECTION_ID> reference_connection_data;

    struct REFERENCE_IRREGULAR_DATA_HELPER
    {
        BLOCK_ID b=BLOCK_ID(-7);
        DOF_PAIRS irreg_pairs[2][2]; // dof[to][from]; 0=regular, 1=irregular
        BLOCK_MATRIX<T> M;
        int mat_id=-7;
    };

    struct REFERENCE_IRREGULAR_DATA
    {
        IRREG_ID ic_id;
        ARRAY<REFERENCE_IRREGULAR_DATA_HELPER> pairs;
    };

    ARRAY<REFERENCE_IRREGULAR_DATA,REFERENCE_IRREGULAR_ID> reference_irregular_data;

    ~COMPONENT_LAYOUT_FEM();
    void Parse_Input(const std::string& pipe_file);

    MATRIX<T,2> Compute_Xform(const TV& dir); // dir is normalized
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint(const JOINT_KEY& key);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint_2(const JOINT_KEY& key);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint_3_Small(const JOINT_KEY& key);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint_3_Average(const JOINT_KEY& key);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint_3(const JOINT_KEY& key);
    CANONICAL_COMPONENT* Make_Canonical_Pipe(const PIPE_KEY& key);
    CANONICAL_BLOCK_ID Make_Canonical_Pipe_Block(const PIPE_KEY& key);
    CANONICAL_COMPONENT* Make_Canonical_Pipe_Change(const PIPE_CHANGE_KEY& key);
    CANONICAL_BLOCK_ID Make_Canonical_Change_Block(const PIPE_CHANGE_KEY& key);
    void Compute();
    void Update_Masters();

    typedef TRIPLE<CANONICAL_BLOCK_ID,INTERVAL<int>,INTERVAL<int> > BC_KEY;
    std::map<PIPE_KEY,BC_KEY> canonical_bc_blocks[2];
    BC_KEY Make_BC_Block(const PIPE_KEY& key,bool is_v);
    
    struct VERTEX_DATA
    {
        TV X;
        BLOCK_CONNECTION con;
        IRREG_ID ic=IRREG_ID(-7); // if con.id<0, vd holds irregular connection. ic is an index in irregular_connections.
    };

    void Emit_Component_Blocks(const CANONICAL_COMPONENT* cc,const XFORM<TV>& xf,ARRAY<VERTEX_DATA>& vd);
    void Set_Connector(VERTEX_DATA& vd,BLOCK_ID id,CON_ID con_id);

    // return: flags indicating which connections interact
    // pair: block + master mask
    HASHTABLE<PAIR<CANONICAL_BLOCK_ID,int>,int> separates_dofs;
    int Separates_Dofs(BLOCK_ID b);

    HASHTABLE<std::tuple<CANONICAL_BLOCK_ID,CON_ID,CANONICAL_BLOCK_ID,CON_ID>,PAIR<CANONICAL_BLOCK_ID,ARRAY<int> > > merge_canonical_blocks;
    void Merge_Blocks(BLOCK_ID id,CON_ID con_id);
    PAIR<CANONICAL_BLOCK_ID,ARRAY<int> >*
        Merge_Canonical_Blocks(CANONICAL_BLOCK_ID id0,CON_ID con_id0,
            XFORM<TV> xf0,CANONICAL_BLOCK_ID id1,CON_ID con_id1,XFORM<TV> xf1);
    int Approx_Dof_Count(BLOCK_ID b);
    void Merge_Blocks();
    void Compute_Matrix_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Fill_Canonical_Block_Matrix(BLOCK_MATRIX<T>& mat,const CANONICAL_BLOCK& cb);
    void Fill_Block_Matrix(REFERENCE_BLOCK_DATA& rd);
    void Fill_Connection_Matrix(REFERENCE_CONNECTION_DATA& cd);
    void Fill_Irregular_Connection_Matrix(REFERENCE_IRREGULAR_DATA& ri);
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
    void Copy_Matrix_Data(BLOCK_MATRIX<T>& A,BLOCK_ID b,
        const DOF_PAIRS& dpa,const DOF_PAIRS& dpb,BLOCK_ID ar,BLOCK_ID ac) const;
    void Copy_Vector_Data(const BLOCK_VECTOR<T>& B,BLOCK_ID b,const DOF_PAIRS& dp,BLOCK_ID a);
    void Init_Block_Matrix(BLOCK_MATRIX<T>& M,BLOCK_ID a,BLOCK_ID b) const;
    void Init_Block_Vector(BLOCK_VECTOR<T>& M,BLOCK_ID b) const;
    void Init_Block_Vector(BLOCK_VECTOR<T>& M,const CANONICAL_BLOCK& cb) const;
    void Times_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& v,const BLOCK_VECTOR<T>& u) const;
    void Times_P_U(BLOCK_ID b,BLOCK_VECTOR<T>& w,const ARRAY<T>& div_v,const ARRAY<T>& div_e) const;
    void Times_Line_Integral_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& w,const BLOCK_VECTOR<T>& u) const;
    void Apply_To_RHS(BLOCK_ID b,const BLOCK_VECTOR<T>& w);
    RANGE<TV> Compute_Bounding_Box() const;
    void Eliminate_Irregular_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Eliminate_Non_Seperators(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Eliminate_Strip(CACHED_ELIMINATION_MATRIX<T>& cem,const ARRAY<BLOCK_ID>& a);
    void Eliminate_Simple(CACHED_ELIMINATION_MATRIX<T>& cem,BLOCK_ID first,CON_ID con_id_source);
    void Visualize_Block_State(BLOCK_ID b) const;
    void Transform_Solution(const CACHED_ELIMINATION_MATRIX<T>& cem);
    void Visualize_Solution(BLOCK_ID b) const;
    void Dump_World_Space_System() const;
    void Transform_To_World_Space(BLOCK_MATRIX<T>& M,const BLOCK_MATRIX<T>& B,BLOCK_ID a,BLOCK_ID b) const;

  private:
    std::tuple<TV,T,T> Elbow_Pit(T angle,T width) const;
    TV Elbow_Pit_Oriented(T angle,T width) const;
    VECTOR<TV,2> Extrude(const TV& v0,const TV& v1,const TV& n) const;
    PAIR<ARRAY<TV>,ARRAY<TV> > Arc(const TV& c,T angle,T len_arm,T ext0,T ext1) const;
    ARRAY<TV> Polyline(const ARRAY<TV>& points,T dx) const;
    void Joint_Connection(int offset,BLOCK_MESHING_ITERATOR<TV>& it,CANONICAL_BLOCK& cb,
        IRREGULAR_CONNECTION& ic0,IRREGULAR_CONNECTION& ic1,ARRAY<BLOCK_CONNECTION,CON_ID>& con,CON_ID prev) const;
};

}
#endif
