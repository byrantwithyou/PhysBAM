//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPONENT_LAYOUT_FEM__
#define __COMPONENT_LAYOUT_FEM__
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <map>
#include "COMMON.h"

namespace PhysBAM{

template<class TV> struct COMPONENT_LAYOUT_FEM;

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

    static constexpr T comp_tol=(T)1e-10;
    T target_length;

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
    
    // ranges of shared dof indices owned by this component and used from
    // neigboring component.
    struct CROSS_SECTION
    {
        INTERVAL<int> owned,used;
        int master_index;
        CROSS_SECTION_TYPE_ID type;
    };

    struct CANONICAL_BLOCK
    {
        ARRAY<CROSS_SECTION> cross_sections;
        ARRAY<TV> X;
        ARRAY<IV3> E;
    };
    ARRAY<CANONICAL_BLOCK,CANONICAL_BLOCK_ID> canonical_blocks;

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
        int num_dofs;
    };

    // regular is master
    struct IRREGULAR_CONNECTION
    {
        BLOCK_ID regular;
        int con_id;
        ARRAY<PAIR<BLOCK_ID,int> > edge_on; // one for each dof on cross section, order: owned,master_index,used
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

    void Parse_Input(const std::string& pipe_file);

    XFORM_ID Compute_Xform(const TV& dir); // dir is normalized
    XFORM Compose_Xform(const XFORM& a,const XFORM& b);
    PAIR<CANONICAL_COMPONENT*,ARRAY<T> > Make_Canonical_Joint(const JOINT_KEY& key);
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
    };

    void Emit_Component_Blocks(CANONICAL_COMPONENT* cc,const XFORM& xf,ARRAY<VERTEX_DATA>& vd);
    void Set_Cornector(VERTEX_DATA& vd,BLOCK_ID id,int con_id);

    // return: flags indicating which connections interact
    // pair: block + master mask
    HASHTABLE<PAIR<CANONICAL_BLOCK_ID,int>,int> separates_dofs;
    int Separates_Dofs(BLOCK_ID b);

    HASHTABLE<std::tuple<CANONICAL_BLOCK_ID,int,CANONICAL_BLOCK_ID,int>,CANONICAL_BLOCK_ID> merge_canonical_blocks;
    void Merge_Blocks(BLOCK_ID id,int con_id);
    CANONICAL_BLOCK_ID Merge_Canonical_Blocks(CANONICAL_BLOCK_ID id0,int con_id0,
        XFORM xf0,CANONICAL_BLOCK_ID id1,int con_id1,XFORM xf1);
    int Approx_Dof_Count(BLOCK_ID b);
    void Merge_Blocks();
};

}
#endif
