//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_LAYOUT_FEM_EXTRUDED__
#define __FLUID_LAYOUT_FEM_EXTRUDED__
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include "COMMON.h"
#include "PARSE_DATA_FEM.h"

namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;
template<class TV> struct FLUID_LAYOUT_FEM;

template<typename T>
struct FLUID_LAYOUT_FEM<VECTOR<T,3> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<T,3> TV3;
    static const TV3 extrude_dir;
    FLUID_LAYOUT_FEM<TV> fl2;
    TETRAHEDRALIZED_VOLUME<T>& vol_hidden;
    ARRAY<typename FLUID_LAYOUT_FEM<TV>::BLOCK_DATA,BLOCK_ID>& blocks=fl2.blocks;
    ARRAY<typename FLUID_LAYOUT_FEM<TV>::PIPE_DATA,PIPE_ID>& pipes=fl2.pipes;
    ARRAY<BLOCK_ID,TRIANGLE_ID> elem_blocks;
    ARRAY<DOF_ID,PARTICLE_ID> pressure_dofs,vel_node_dofs;
    ARRAY<DOF_ID,EDGE_ID> vel_edge_dofs;
    ARRAY<VECTOR<int,2>,DOF_ID> dof_map;
    DOF_ID num_dofs=DOF_ID(0);

    int layers=-1;
    BC_ID wall_bc=BC_ID(-1);

    FLUID_LAYOUT_FEM();
    ~FLUID_LAYOUT_FEM();
    void Dump_Input(const PARSE_DATA_FEM<TV,TV3>& pd) const;
    void Dump_Mesh() const;
    void Dump_Layout() const;
    void Dump_Node_Blocks() const;
    void Dump_Edge_Blocks() const;
    void Dump_Dofs() const;
    void Dump_Blocks() const;
    void Print_Statistics() const;
    void Extrude(int layers,T dz);
    void Update_BC_Planar(const PARSE_DATA_FEM<TV,TV3>& pd);
    void Assign_Edge_Blocks_Planar();
    void Allocate_Dofs(const PARSE_DATA_FEM<TV,TV3>& pd);
    void Compute(const PARSE_DATA_FEM<TV,TV3>& pd);

    BC_ID Node_BC(PARTICLE_ID p) const
    {
        if(Value(p)<fl2.area_hidden.particles.number || Value(p)>=layers*fl2.area_hidden.particles.number)
            return wall_bc;
        BC_ID r(-1);
        PARTICLE_ID q(Value(p)%fl2.area_hidden.particles.number);
        const BC_ID* pbc=fl2.particle_bc_map.Get_Pointer(q);
        if(pbc) r=*pbc;
        return r;
    }
    BC_ID Edge_BC(EDGE_ID e) const
    {
        VECTOR<PARTICLE_ID,2> p=Edge(e);
        int m=fl2.area_hidden.particles.number;
        if((Value(p(0))<m && Value(p(1))<m) || (Value(p(0))>=layers*m && Value(p(1))>=layers*m))
            return wall_bc;
        BC_ID r(-1);
        const BC_ID* ebc=0;
        VECTOR<PARTICLE_ID,2> q={PARTICLE_ID(Value(p(0))%m),PARTICLE_ID(Value(p(1))%m)};
        if(q(0)==q(1))
            ebc=fl2.particle_bc_map.Get_Pointer(q(0));
        else
            ebc=fl2.bc_map.Get_Pointer(q.Sorted());
        if(ebc) r=*ebc;
        return r;
    }
    const TV3& X(PARTICLE_ID p) const {return vol_hidden.particles.X(Value(p));}
    PARTICLE_ID Number_Particles() const {return PARTICLE_ID(vol_hidden.particles.number);}
    TRIANGLE_ID Number_Tetrahedrons() const {return TRIANGLE_ID(vol_hidden.mesh.elements.m);}
    EDGE_ID Number_Edges() const {return EDGE_ID(vol_hidden.mesh.segment_mesh->elements.m);}
    VECTOR<PARTICLE_ID,2> Edge(EDGE_ID e) const
    {
        VECTOR<int,2> v=vol_hidden.mesh.segment_mesh->elements(Value(e));
        return {PARTICLE_ID(v(0)),PARTICLE_ID(v(1))};
    }
    VECTOR<int,4> Tetrahedron(TRIANGLE_ID t) const
    {
        auto a=vol_hidden.mesh.elements(Value(t));
        return {a(0),a(1),a(2),a(3)};
    }
    BLOCK_ID Node_Block(PARTICLE_ID p) const
    {
        PARTICLE_ID q(Value(p)%fl2.area_hidden.particles.number);
        return fl2.node_blocks(q);
    }
    BLOCK_ID Edge_Block(PARTICLE_ID p0,PARTICLE_ID p1) const
    {
        PARTICLE_ID q0(Value(p0)%fl2.area_hidden.particles.number),
            q1(Value(p1)%fl2.area_hidden.particles.number);
        if(q0==q1) return fl2.node_blocks(q0);
        else return fl2.edge_blocks(EDGE_ID(fl2.area_hidden.mesh.segment_mesh->Segment(Value(q0),Value(q1))));
    }
};

}
#endif
