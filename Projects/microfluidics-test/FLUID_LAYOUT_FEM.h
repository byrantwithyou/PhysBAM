//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_LAYOUT_FEM__
#define __FLUID_LAYOUT_FEM__
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include "COMMON.h"
#include "PARSE_DATA_FEM.h"

namespace PhysBAM{

template<class T> class TRIANGULATED_AREA;
template<class TV> struct FLUID_LAYOUT_FEM;

template<class T>
struct FLUID_LAYOUT_FEM<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;

    // (vert index,pipe index) -> array of particle id
    typedef HASHTABLE<PAIR<VERTEX_ID,PIPE_ID>,ARRAY<PARTICLE_ID> > CONNECTION;

    struct ELEMENT_DATA
    {
        BLOCK_ID block_id;
        int greed;
    };

    struct BLOCK_DATA
    {
        PIPE_ID pipe_id;
        int num_dofs;
    };

    struct PIPE_DATA
    {
        TV u,v; // u = x1-x0, v=x2-x0
        TRIANGLE_ID first_element;
    };

    TRIANGULATED_AREA<T>& area_hidden;

    ARRAY<BLOCK_DATA,BLOCK_ID> blocks;
    ARRAY<ELEMENT_DATA,TRIANGLE_ID> elem_data; // per element
    HASHTABLE<VECTOR<PARTICLE_ID,2>,BC_ID> bc_map;
    HASHTABLE<PARTICLE_ID,BC_ID> particle_bc_map;
    ARRAY<PIPE_DATA,PIPE_ID> pipes; // per pipe
    
    ARRAY<BLOCK_ID,PARTICLE_ID> node_blocks;
    ARRAY<DOF_ID,PARTICLE_ID> pressure_dofs,vel_node_dofs;
    ARRAY<BLOCK_ID,EDGE_ID> edge_blocks;
    ARRAY<DOF_ID,EDGE_ID> vel_edge_dofs;
    ARRAY<VECTOR<int,2>,DOF_ID> dof_map;
    DOF_ID num_dofs;
    
    FLUID_LAYOUT_FEM();
    ~FLUID_LAYOUT_FEM();
    template<class TW>
    void Compute(const PARSE_DATA_FEM<TV,TW>& pd);
    template<class TW>
    void Allocate_Dofs(const PARSE_DATA_FEM<TV,TW>& pd);
    void Print_Statistics() const;
    template<class TW>
    void Generate_End(VERTEX_ID i,PIPE_ID pipe,const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con);
    template<class TW>
    void Generate_Joint(VERTEX_ID i,const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con);
    template<class TW>
    void Generate_2_Joint(VERTEX_ID i,const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con,JOINT_TYPE jt);
    template<class TW>
    void Generate_3_Joint(VERTEX_ID i,const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con);
    template<class TW>
    void Generate_3_Joint_SmallMin(VERTEX_ID i,const VECTOR<VERTEX_ID,3>& ends,const VECTOR<PIPE_ID,3>& pipes,
        const VECTOR<T,3>& angles,const VECTOR<VECTOR<TV,3>,3>& tri,
        const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con);
    template<class TW>
    void Generate_3_Joint_LargeMin(VERTEX_ID i,const VECTOR<VERTEX_ID,3>& ends,const VECTOR<PIPE_ID,3>& pipes,
        const VECTOR<T,3>& angles,const VECTOR<VECTOR<TV,3>,3>& tri,
        const PARSE_DATA_FEM<TV,TW>& pd,CONNECTION& con);
    template<class TW>
    void Generate_Pipe(PIPE_ID pipe,const PARSE_DATA_FEM<TV,TW>& pd,const CONNECTION& con);

    void Dump_Mesh() const;
    void Dump_Layout() const;
    template<class TW>
    void Dump_Input(const PARSE_DATA_FEM<TV,TW>& pd) const;
    void Dump_Edge_Blocks() const;
    void Dump_Node_Blocks() const;
    void Dump_Dofs() const;

    void Mark_BC(const ARRAY<PARTICLE_ID>& pindices,BC_ID bc_id);
    // return (center, normalized start vec, normalied end vec)
    VECTOR<TV,3> Wedge(const TV& joint,const TV& p0,const TV& p1,int half_width,T unit_length) const;
    ARRAY<PARTICLE_ID> Sample_Interpolated(T s,const ARRAY<PARTICLE_ID>& side0,
        const ARRAY<PARTICLE_ID>& side1,T unit_length);
    void Merge_Interpolated(const ARRAY<PARTICLE_ID>& left,const ARRAY<PARTICLE_ID>& right);
    ARRAY<PARTICLE_ID> Polyline(const ARRAY<TV>& points,T unit_length);
    PAIR<ARRAY<PARTICLE_ID>,ARRAY<PARTICLE_ID> > Arc(const TV& c,
        const TV& p0,const TV& p1,int half_width,T unit_length,const TV& dir0,
        T n0,const TV& dir1,T n1);
    PAIR<ARRAY<PARTICLE_ID>,ARRAY<PARTICLE_ID> > Corner(const TV& c,const TV& joint,
        const TV& p0,const TV& p1,T unit_length,const TV& dir0,T n0,const TV& dir1,T n1);
    void Weld(int n,const ARRAY<PARTICLE_ID>& side0,const ARRAY<PARTICLE_ID>& side1,
        T unit_length,ARRAY<PARTICLE_ID>& f0,ARRAY<PARTICLE_ID>& f1);

    const TV& X(PARTICLE_ID p) const {return area_hidden.particles.X(Value(p));}
    TV& X(PARTICLE_ID p) {return area_hidden.particles.X(Value(p));}
    PARTICLE_ID Add_Particle()
    {
        node_blocks.Append(BLOCK_ID(-1));
        return PARTICLE_ID(area_hidden.particles.Add_Element());
    }
    PARTICLE_ID Add_Particles(int n)
    {
        node_blocks.Resize(node_blocks.m+n,use_init,BLOCK_ID(-1));
        return PARTICLE_ID(area_hidden.particles.Add_Elements(n));
    }
    PARTICLE_ID Number_Particles() const {return PARTICLE_ID(area_hidden.particles.number);}
    TRIANGLE_ID Number_Triangles() const {return TRIANGLE_ID(area_hidden.mesh.elements.m);}
    EDGE_ID Number_Edges() const {return EDGE_ID(area_hidden.mesh.segment_mesh->elements.m);}

    void Append_Triangle(const VECTOR<PARTICLE_ID,3>& e)
    {area_hidden.mesh.elements.Append({Value(e.x),Value(e.y),Value(e.z)});}

    ARRAY_VIEW<const TRIANGLE_ID> Edge_Triangles(EDGE_ID t) const
    {
        const auto& a=(*area_hidden.mesh.edge_triangles)(Value(t));
        return {reinterpret_cast<const TRIANGLE_ID*>(a.Get_Array_Pointer()),a.m};
    }
    ARRAY_VIEW<const TRIANGLE_ID> Incident_Triangles(PARTICLE_ID t) const
    {
        const auto& a=(*area_hidden.mesh.incident_elements)(Value(t));
        return {reinterpret_cast<const TRIANGLE_ID*>(a.Get_Array_Pointer()),a.m};
    }
    VECTOR<PARTICLE_ID,2> Edge(EDGE_ID e) const
    {
        auto a=area_hidden.mesh.segment_mesh->elements(Value(e));
        return {PARTICLE_ID(a(0)),PARTICLE_ID(a(1))};
    }
    VECTOR<PARTICLE_ID,3> Triangle(TRIANGLE_ID t) const
    {
        auto a=area_hidden.mesh.elements(Value(t));
        return {PARTICLE_ID(a(0)),PARTICLE_ID(a(1)),PARTICLE_ID(a(2))};
    }
    VECTOR<EDGE_ID,3> Triangle_Edges(TRIANGLE_ID t) const
    {
        auto& a=(*area_hidden.mesh.element_edges)(Value(t));
        return {EDGE_ID(a(0)),EDGE_ID(a(1)),EDGE_ID(a(2))};
    }
    T Area(TRIANGLE_ID t) const
    {
        return area_hidden.Area(Value(t));
    }
};

}
#endif
