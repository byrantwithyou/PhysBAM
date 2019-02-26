//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPONENT_LAYOUT_FEM__
#define __COMPONENT_LAYOUT_FEM__
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_VECTOR.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <map>
#include "BLOCK_MATRIX.h"
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

    enum class BC_TYPE {dirichlet_v,traction,analytic};
    struct BOUNDARY_CONDITION
    {
        CROSS_SECTION_TYPE_ID type;
        BC_TYPE bc_type;
        TV X,normal;
        T flowrate;
        TV traction;
    };
    ANALYTIC_VECTOR<TV>* force=0,*analytic_velocity=0;
    ANALYTIC_SCALAR<TV>* analytic_pressure=0;
    ARRAY<BOUNDARY_CONDITION,BC_ID> boundary_conditions;

    TV Velocity(BC_ID id,const TV& X) const
    {
        if(analytic_velocity) return analytic_velocity->v(X,0);
        else
        {
            const BOUNDARY_CONDITION& bc=boundary_conditions(id);
            const auto& cst=cross_section_types(bc.type);
            T a=cst.width/2;
            T r=(X-bc.X).Magnitude();
            T v=-3.0/4*bc.flowrate/(a*a*a)*(r*r-a*a);
            return v*bc.normal;
        }
    }

    TV Traction(BC_ID id,const TV& X) const
    {
        const BOUNDARY_CONDITION& bc=boundary_conditions(id);
        if(analytic_velocity && analytic_pressure)
        {
            SYMMETRIC_MATRIX<T,TV::m> stress=analytic_velocity->dX(X,0).Twice_Symmetric_Part()*mu;
            stress-=analytic_pressure->f(X,0);
            return stress*bc.normal;
        }
        else return bc.traction;
    }

    TV Force(const TV& X) const
    {
        if(analytic_velocity && analytic_pressure)
        {
            SYMMETRIC_TENSOR<T,0,TV::m> ddU=analytic_velocity->ddX(X,0);
            TV f=analytic_pressure->dX(X,0);
            f-=mu*(Contract<1,2>(ddU)+Contract<0,2>(ddU));
            return f;
        }
        else return force?force->v(X,0):TV();
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

    typedef std::tuple<CANONICAL_BLOCK_ID,int,CANONICAL_BLOCK_ID,int> REGULAR_CON_KEY;
    HASHTABLE<REGULAR_CON_KEY,int> regular_connection_matrix_blocks;
    
    // first is master, second is slave
    struct BLOCK_CONNECTION
    {
        BLOCK_ID id;
        int con_id; // if irregular, ~con_id is index into irregular_connections
        bool master;
    };

    struct BLOCK
    {
        CANONICAL_BLOCK_ID block;
        XFORM xform;
        ARRAY<BLOCK_CONNECTION> connections;
        ARRAY<int> edge_on; // for edge-on (index in irregular_connections)
    };

    // regular is master
    struct IRREGULAR_CONNECTION
    {
        BLOCK_ID regular;
        int con_id;
        // one for each dof on cross section, starting from owned side of cross section
        ARRAY<PAIR<BLOCK_ID,int> > edge_on_v,edge_on_e;
        int ref_ic;
        int block_data;
    };

    // neighbor block i is given index ~i and con_id=-1.
    struct CANONICAL_COMPONENT
    {
        ARRAY<BLOCK,BLOCK_ID> blocks;
        ARRAY<BLOCK_CONNECTION> connections;
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

    struct REFERENCE_BLOCK_DATA
    {
        int num_dofs_v,num_dofs_e;
        ARRAY<int> dof_map_v,dof_map_e;
    };

    ARRAY<REFERENCE_BLOCK_DATA> reference_block_data;

    struct IRREGULAR_REFERENCE_BLOCK_DATA
    {
        ARRAY<int> add_block;
    };

    ARRAY<IRREGULAR_REFERENCE_BLOCK_DATA> irregular_reference_block_data;

    ~COMPONENT_LAYOUT_FEM();
    void Parse_Input(const std::string& pipe_file);

    XFORM_ID Compute_Xform(const TV& dir); // dir is normalized
    XFORM Compose_Xform(const XFORM& a,const XFORM& b);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint(const JOINT_KEY& key);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint_2(const JOINT_KEY& key);
    CANONICAL_COMPONENT* Make_Canonical_Pipe(const PIPE_KEY& key);
    CANONICAL_BLOCK_ID Make_Canonical_Pipe_Block(const PIPE_KEY& key);
    CANONICAL_COMPONENT* Make_Canonical_Pipe_Change(const PIPE_CHANGE_KEY& key);
    CANONICAL_BLOCK_ID Make_Canonical_Change_Block(const PIPE_CHANGE_KEY& key);
    CROSS_SECTION_TYPE_ID Get_Cross_Section_ID(const CROSS_SECTION_TYPE& cs);
    void Compute();
    void Update_Masters();

    struct VERTEX_DATA
    {
        TV X;
        IRREGULAR_CONNECTION con;
        BC_ID bc_id=BC_ID(-1);
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
        const ARRAY<IV>& va,const ARRAY<IV>& ea,
        const ARRAY<IV>& vb,const ARRAY<IV>& eb,
        BLOCK_ID ar,BLOCK_ID ac) const;
    void Init_Block_Matrix(BLOCK_MATRIX<T>& M,BLOCK_ID a,BLOCK_ID b) const;
    void Compute_Reference_Irregular_Connections();
  private:
    std::tuple<TV,T,T> Vertex(T angle,T width) const;
    PAIR<ARRAY<TV>,ARRAY<TV> > Arc(const TV& c,T angle,T len_arm,T ext0,T ext1) const;
    ARRAY<IV3> Merge_Interpolated(const ARRAY<TV>& X,int n0,int n1) const;
    ARRAY<TV> Interpolated(T s,const ARRAY<TV>& side0,const ARRAY<TV>& side1) const;
    ARRAY<std::tuple<ARRAY<TV>,ARRAY<IV3>,int> > Fill(int nseg,const ARRAY<TV>& inner,const ARRAY<TV>& outer) const;
};

}
#endif
