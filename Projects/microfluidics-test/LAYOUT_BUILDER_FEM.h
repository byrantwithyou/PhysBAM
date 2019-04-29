//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LAYOUT_BUILDER_FEM__
#define __LAYOUT_BUILDER_FEM__
#include "COMPONENT_BC.h"
#include "COMPONENT_CHANGE.h"
#include "COMPONENT_JOINT.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include "COMPONENT_PIPE.h"

namespace PhysBAM{
template<class T>
struct LAYOUT_BUILDER_FEM
{
    typedef VECTOR<T,2> TV;
    PHYSBAM_DECLARE_ELEMENT_ID(CS_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
    PHYSBAM_DECLARE_ELEMENT_ID(VERT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
    PHYSBAM_DECLARE_ELEMENT_ID(CONNECTOR_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
    PHYSBAM_DECLARE_ELEMENT_ID(BC_U_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
    PHYSBAM_DECLARE_ELEMENT_ID(BC_T_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
    enum class TAG {SET_DEPTH,SET_LEN,DECL_CS,DECL_VERT,SET_BC_U,SET_BC_T,DECL_JT,DECL_PIPE,CS,VERT,CONNECTOR};
    struct VERTEX_DATA
    {
        TV X;
        BLOCK_CONNECTION con;
    };

    COMPONENT_LAYOUT_FEM<T>& cl;
    ARRAY<PAIR<TAG,int> > commands;
    ARRAY<PAIR<int,T>,CS_ID> cross_sections;
    ARRAY<TV,VERT_ID> verts;
    ARRAY<VECTOR<T,2> > var_size_pipes; // (distance, length)
    CONNECTOR_ID last_cid=CONNECTOR_ID(0);
    HASHTABLE<CONNECTOR_ID,VERTEX_DATA> connectors;

    COMPONENT_PIPE<T> comp_pipe;
    COMPONENT_CHANGE<T> comp_change;
    COMPONENT_BC<T> comp_bc;
    COMPONENT_JOINT<T> comp_joint;

    LAYOUT_BUILDER_FEM(COMPONENT_LAYOUT_FEM<T>& cl);
    void Emit_Component_Blocks(const CANONICAL_COMPONENT<T>* cc,const XFORM<TV>& xf,ARRAY<VERTEX_DATA>& vd);
    void Set_Connector(VERTEX_DATA& vd,BLOCK_ID id,CON_ID con_id);
    MATRIX<T,2> Compute_Xform(const TV& dir); // dir is normalized
    // l <target-length>
    // z <depth> <volumetric-layers>
    // c <cross-section-name> <num-elements> <width>
    // v <vertex-name> <vertex-location-2d>
    // j <cross-section-name> <num-pipes> <origin-vertex> [<vertex-name> <connection-name>]*
    // p <cross-section-name> <connection-name> <connection-name>
    // g <cross-section-name> <cross-section-name> <vertex-name> <vertex-name> <distance> <length> <connection-name> <connection-name>
    // u <cross-section-name> <origin-vertex> <vertex-name> <connection-name> <flow-rate>
    // t <cross-section-name> <origin-vertex> <vertex-name> <connection-name> <traction-2d>
    void From_File(const std::string& file);
    std::string To_String() const;

    void Set_Target_Length(T l);
    void Set_Depth(T z,int m);
    CS_ID Cross_Section(int d,T w);
    VERT_ID Vertex(const TV& X);
    PAIR<CONNECTOR_ID,BC_U_ID> Set_BC(CS_ID cs,VERT_ID from,VERT_ID to,T flow_rate);
    PAIR<CONNECTOR_ID,BC_T_ID> Set_BC(CS_ID cs,VERT_ID from,VERT_ID to,const TV& traction);
    ARRAY<CONNECTOR_ID> Joint(CS_ID cs,int n,VERT_ID o,const ARRAY<VERT_ID>& arms);
    void Pipe(CS_ID cs,CONNECTOR_ID a,CONNECTOR_ID b);
    VECTOR<CONNECTOR_ID,2> Pipe(CS_ID cs0,CS_ID cs1,VERT_ID v0,VERT_ID v1,T offset,T length);
};
}
#endif
