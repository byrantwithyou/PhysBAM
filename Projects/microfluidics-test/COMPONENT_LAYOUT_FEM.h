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

namespace PhysBAM{

template<class TV> struct COMPONENT_LAYOUT_FEM;
template<class T> struct CACHED_ELIMINATION_MATRIX;

PHYSBAM_DECLARE_ELEMENT_ID(XFORM_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(COMPONENT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(SEPARATOR_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(CANONICAL_BLOCK_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(CROSS_SECTION_TYPE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);

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
    
    ARRAY<MATRIX<T,TV::m>,XFORM_ID> xforms;
    HASHTABLE<PAIR<XFORM_ID,XFORM_ID>,XFORM_ID> xform_comp_table;
    
    struct XFORMS_LOOKUP_CMP
    {
        bool operator()(const TV& a,const TV& b)
        {
            for(int i=0;i<TV::m;i++)
            {
                if(a(i)<b(i)) return true;
                if(b(i)<a(i)) return false;
            }
            return false;
        }
    };

    std::map<TV,XFORM_ID,XFORMS_LOOKUP_CMP> xforms_lookup;

    // transform (rotation + uniform scale + shift); y=M*x+b;
    struct XFORM
    {
        XFORM_ID id; // for matrix part
        TV b;
    };

    struct CROSS_SECTION_TYPE
    {
        int num_dofs;
        T width;
        int Cmp(const CROSS_SECTION_TYPE& c) const
        {
            if(num_dofs<c.num_dofs) return -1;
            if(c.num_dofs<num_dofs) return 1;
            if(width<c.width-comp_tol) return -1;
            if(c.width<width-comp_tol) return 1;
            return 0;
        }

        bool operator<(const CROSS_SECTION_TYPE& c) const
        {return Cmp(c)<0;}
    };

    ARRAY<CROSS_SECTION_TYPE,CROSS_SECTION_TYPE_ID> cross_section_types;
    std::map<CROSS_SECTION_TYPE,CROSS_SECTION_TYPE_ID> cross_section_type_lookup;
    
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
        ARRAY<CROSS_SECTION> cross_sections;
        ARRAY<TV> X;
        ARRAY<IV3> E;
        ARRAY<IV> S;
        ARRAY<int> bc_v,bc_e;
    };
    ARRAY<CANONICAL_BLOCK,CANONICAL_BLOCK_ID> canonical_blocks;

    ARRAY<BLOCK_MATRIX<T>,CANONICAL_BLOCK_ID> canonical_block_matrices;

    ARRAY<BLOCK_MATRIX<T> > matrix_block_list;

    ARRAY<BLOCK_VECTOR<T>,BLOCK_ID> rhs_block_list;
    
    typedef std::tuple<CANONICAL_BLOCK_ID,int,CANONICAL_BLOCK_ID,int> REGULAR_CON_KEY;
    HASHTABLE<REGULAR_CON_KEY,int> regular_connection_matrix_blocks;
    
    // first is master, second is slave
    struct BLOCK_CONNECTION
    {
        BLOCK_ID id=BLOCK_ID(-7);;
        int con_id; // if irregular, ~con_id is index into irregular_connections
        bool master;
    };

    struct BLOCK
    {
        CANONICAL_BLOCK_ID block;
        XFORM xform;
        ARRAY<BLOCK_CONNECTION> connections;
        ARRAY<int> edge_on; // for edge-on (index in irregular_connections)
        int flags=0; // 1=separator, 2=separator-eligible
    };

    // regular is master
    struct IRREGULAR_CONNECTION
    {
        BLOCK_ID regular=BLOCK_ID(-7);;
        int con_id;
        // one for each dof on cross section, starting from owned side of cross section
        ARRAY<PAIR<BLOCK_ID,int> > edge_on_v,edge_on_e;
        int ref_ic=-7;
        int block_data=-7;
    };

    // neighbor block i is given index ~i and con_id=-1.
    struct CANONICAL_COMPONENT
    {
        ARRAY<BLOCK,BLOCK_ID> blocks;
        ARRAY<IRREGULAR_CONNECTION> irregular_connections;
    };

    ARRAY<BLOCK,BLOCK_ID> blocks;
    ARRAY<IRREGULAR_CONNECTION> irregular_connections;
    
    struct PIPE_KEY
    {
        CROSS_SECTION_TYPE_ID type;
        T length;
        bool operator<(const PIPE_KEY& p) const
        {
            if(type!=p.type) return type<p.type;
            if(length<p.length-comp_tol) return true;
            return false;
        }
    };
    std::map<PIPE_KEY,CANONICAL_COMPONENT*> canonical_pipes;
    std::map<PIPE_KEY,CANONICAL_BLOCK_ID> canonical_pipe_blocks;

    struct JOINT_KEY
    {
        CROSS_SECTION_TYPE_ID type;
        ARRAY<T> angles;
        bool operator<(const JOINT_KEY& p) const
        {
            if(type!=p.type) return type<p.type;
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
        CROSS_SECTION_TYPE_ID type[2];
        T length;

        bool operator<(const PIPE_CHANGE_KEY& p) const
        {
            for(int i=0;i<2;i++)
                if(type[i]!=p.type[i])
                    return type[i]<p.type[i];
            if(length<p.length-comp_tol) return true;
            return false;
        }
    };
    std::map<PIPE_CHANGE_KEY,CANONICAL_COMPONENT*> canonical_changes;
    std::map<PIPE_CHANGE_KEY,CANONICAL_BLOCK_ID> canonical_change_blocks;

    ARRAY<PAIR<BLOCK_ID,int>,BLOCK_ID> reference_block; // block, index

    struct DOF_PAIRS
    {
        ARRAY<IV> v,e,p;
    };
    
    struct REFERENCE_BLOCK_DATA
    {
        int num_dofs_v=-7,num_dofs_e=-7,num_dofs_p=-7;
        ARRAY<int> dof_map_v,dof_map_e,dof_map_p;
        DOF_PAIRS pairs;
        ARRAY<DOF_PAIRS> regular_pairs; // this block is the source, the connection is the destination
    };

    ARRAY<REFERENCE_BLOCK_DATA> reference_block_data;

    struct IRREGULAR_REFERENCE_BLOCK_DATA_HELPER
    {
        BLOCK_ID b=BLOCK_ID(-7);;
        int mat_id=-7;
        DOF_PAIRS dof[2][2]; // dof[to][from]; 0=regular, 1=irregular
    };

    struct IRREGULAR_REFERENCE_BLOCK_DATA
    {
        ARRAY<IRREGULAR_REFERENCE_BLOCK_DATA_HELPER> pairs;
    };

    ARRAY<IRREGULAR_REFERENCE_BLOCK_DATA> irregular_reference_block_data;

    ~COMPONENT_LAYOUT_FEM();
    void Parse_Input(const std::string& pipe_file);

    XFORM_ID Compute_Xform(const TV& dir); // dir is normalized
    XFORM Compose_Xform(const XFORM& a,const XFORM& b);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint(const JOINT_KEY& key);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint_2(const JOINT_KEY& key);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint_3_Small(const JOINT_KEY& key);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint_3_Average(const JOINT_KEY& key);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint_3(const JOINT_KEY& key);
    CANONICAL_COMPONENT* Make_Canonical_Pipe(const PIPE_KEY& key);
    CANONICAL_BLOCK_ID Make_Canonical_Pipe_Block(const PIPE_KEY& key);
    CANONICAL_COMPONENT* Make_Canonical_Pipe_Change(const PIPE_CHANGE_KEY& key);
    CANONICAL_BLOCK_ID Make_Canonical_Change_Block(const PIPE_CHANGE_KEY& key);
    CROSS_SECTION_TYPE_ID Get_Cross_Section_ID(const CROSS_SECTION_TYPE& cs);
    void Compute();
    void Update_Masters();

    typedef TRIPLE<CANONICAL_BLOCK_ID,INTERVAL<int>,INTERVAL<int> > BC_KEY;
    std::map<PIPE_KEY,BC_KEY> canonical_bc_blocks[2];
    BC_KEY Make_BC_Block(const PIPE_KEY& key,bool is_v);
    
    struct VERTEX_DATA
    {
        TV X;
        IRREGULAR_CONNECTION con;
    };

    void Emit_Component_Blocks(const CANONICAL_COMPONENT* cc,const XFORM& xf,ARRAY<VERTEX_DATA>& vd);
    void Set_Connector(VERTEX_DATA& vd,BLOCK_ID id,int con_id);

    // return: flags indicating which connections interact
    // pair: block + master mask
    HASHTABLE<PAIR<CANONICAL_BLOCK_ID,int>,int> separates_dofs;
    int Separates_Dofs(BLOCK_ID b);

    HASHTABLE<std::tuple<CANONICAL_BLOCK_ID,int,CANONICAL_BLOCK_ID,int>,PAIR<CANONICAL_BLOCK_ID,ARRAY<int> > > merge_canonical_blocks;
    void Merge_Blocks(BLOCK_ID id,int con_id);
    PAIR<CANONICAL_BLOCK_ID,ARRAY<int> >*
        Merge_Canonical_Blocks(CANONICAL_BLOCK_ID id0,int con_id0,
            XFORM xf0,CANONICAL_BLOCK_ID id1,int con_id1,XFORM xf1);
    int Approx_Dof_Count(BLOCK_ID b);
    void Merge_Blocks();
    void Compute_Matrix_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Fill_Canonical_Block_Matrix(BLOCK_MATRIX<T>& mat,const CANONICAL_BLOCK& cb);
    void Fill_Block_Matrix(BLOCK_ID b);
    int Fill_Connection_Matrix(BLOCK_ID b0,int con_id0,BLOCK_ID b1,int con_id1);
    void Fill_Irregular_Connection_Matrix(IRREGULAR_CONNECTION& ic);
    void Compute_Reference_Blocks();
    int Compute_Connection_Hash(BLOCK_ID b0,int con_id0,BLOCK_ID b1,int con_id1);
    PAIR<int,int> Remap_Owned_Dofs(ARRAY<int>& map_v,ARRAY<int>& map_e,BLOCK_ID b);
    void Compute_Dof_Remapping(BLOCK_ID b);
    void Copy_Matrix_Data(BLOCK_MATRIX<T>& A,BLOCK_ID b,
        const DOF_PAIRS& dpa,const DOF_PAIRS& dpb,BLOCK_ID ar,BLOCK_ID ac) const;
    void Copy_Vector_Data(const BLOCK_VECTOR<T>& B,BLOCK_ID b,const DOF_PAIRS& dp,BLOCK_ID a);
    void Init_Block_Matrix(BLOCK_MATRIX<T>& M,BLOCK_ID a,BLOCK_ID b) const;
    void Init_Block_Vector(BLOCK_VECTOR<T>& M,BLOCK_ID b) const;
    void Init_Block_Vector(BLOCK_VECTOR<T>& M,const CANONICAL_BLOCK& cb) const;
    void Compute_Reference_Irregular_Connections();
    void Times_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& v,const BLOCK_VECTOR<T>& u) const;
    void Times_P_U(BLOCK_ID b,BLOCK_VECTOR<T>& w,const ARRAY<T>& div_v,const ARRAY<T>& div_e) const;
    void Times_Line_Integral_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& w,const BLOCK_VECTOR<T>& u) const;
    void Apply_To_RHS(BLOCK_ID b,const BLOCK_VECTOR<T>& w);
    RANGE<TV> Compute_Bounding_Box() const;
    void Eliminate_Irregular_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Eliminate_Non_Seperators(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Eliminate_Strip(CACHED_ELIMINATION_MATRIX<T>& cem,const ARRAY<BLOCK_ID>& a);
    void Eliminate_Simple(CACHED_ELIMINATION_MATRIX<T>& cem,BLOCK_ID first,int con_id_source);
    void Visualize_Block_State(BLOCK_ID b);

  private:
    std::tuple<TV,T,T> Elbow_Pit(T angle,T width) const;
    TV Elbow_Pit_Oriented(T angle,T width) const;
    VECTOR<TV,2> Extrude(const TV& v0,const TV& v1,const TV& n) const;
    PAIR<ARRAY<TV>,ARRAY<TV> > Arc(const TV& c,T angle,T len_arm,T ext0,T ext1) const;
    ARRAY<TV> Polyline(const ARRAY<TV>& points,T dx) const;
};

}
#endif
