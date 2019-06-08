//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLUTION_FEM__
#define __SOLUTION_FEM__

#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <Geometry/Topology/TETRAHEDRON_MESH.h>
#include <Geometry/Topology/TOPOLOGY_POLICY.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include "BLOCK_VECTOR.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include "MATRIX_CONSTRUCTION_FEM.h"
#include "VISITORS_FEM.h"

namespace PhysBAM{

template<class TV> struct HIERARCHY_POLICY;
template<class T>
struct HIERARCHY_POLICY<VECTOR<T,2> >
{
    typedef TRIANGLE_HIERARCHY_2D<T> TREE;
};
template<class T>
struct HIERARCHY_POLICY<VECTOR<T,3> >
{
    typedef TETRAHEDRON_HIERARCHY<T> TREE;
};

template<class TV>
struct SOLUTION_FEM
{
    typedef typename TV::SCALAR T;
    typedef typename MESH_POLICY<TV::m>::MESH MESH;
    typedef typename HIERARCHY_POLICY<TV>::TREE TREE;

    GEOMETRY_PARTICLES<TV> particles;
    MESH mesh;
    ARRAY<VECTOR<int,(TV::m-1)*3> > element_edges;
    ARRAY<int,BLOCK_ID> first_v,first_e,last_elem;
    ARRAY<PAIR<BLOCK_ID,int> > dof_v,dof_e,dof_p;
    HASHTABLE<PAIR<BLOCK_ID,int>,TV> bc_v,bc_e;
    TREE* tree=0;

    SOLUTION_FEM()
    {
        particles.Store_Velocity(false);
    }

    ~SOLUTION_FEM()
    {
        delete tree;
    }

    void Copy_Mesh_Data(const COMPONENT_LAYOUT_FEM<T>& cl,BLOCK_ID a,BLOCK_ID b,const DOF_PAIRS& dp)
    {
        DOF_LAYOUT<TV> dl0(cl,cl.reference_block_data(cl.blocks(b).ref_id),true);
        DOF_LAYOUT<TV> dl1(cl,cl.reference_block_data(cl.blocks(a).ref_id),false);
        Visit_Dof_Pairs(dl0,dl1,dp,
            [=](int d,int s){dof_v(first_v(a)+s)={b,d};},
            [=](int d,int s){dof_e(first_e(a)+s)={b,d};},
            [=](int d,int s){dof_p(first_v(a)+s)={b,d};});
    }

    void Fill_Global_Dof_Mapping(const COMPONENT_LAYOUT_FEM<T>& cl)
    {
        dof_v.Resize(particles.number,use_init,{BLOCK_ID(-7),-7});
        dof_p.Resize(particles.number,use_init,{BLOCK_ID(-7),-7});
        dof_e.Resize(element_edges.m*(TV::m-1)*3,use_init,{BLOCK_ID(-7),-7});

        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            const auto& bl=cl.blocks(b);
            const auto& rd=cl.reference_block_data(bl.ref_id);

            Copy_Mesh_Data(cl,b,b,rd.pairs);

            for(CON_ID cc(0);cc<bl.connections.m;cc++)
            {
                const auto& c=bl.connections(cc);
                if(c.is_regular)
                    Copy_Mesh_Data(cl,b,c.id,cl.Regular_Connection_Pair(b,cc,false));
                else
                {
                    const auto& ic=cl.irregular_connections(c.irreg_id);
                    const auto& irbd=cl.reference_irregular_data(ic.ref_id);
                    for(int e=0;e<ic.edge_on.m;e++)
                    {
                        if(irbd.mapping(e).y)
                        {
                            const auto& h=irbd.pairs(irbd.mapping(e).x);
                            Copy_Mesh_Data(cl,b,ic.edge_on(e).b,h.irreg_pairs[1]);
                        }
                    }
                }
            }

            for(auto e:bl.edge_on)
            {
                const auto& ic=cl.irregular_connections(e.x);
                const auto& irbd=cl.reference_irregular_data(ic.ref_id);
                if(irbd.mapping(e.y).y)
                {
                    const auto& h=irbd.pairs(irbd.mapping(e.y).x);
                    Copy_Mesh_Data(cl,b,ic.regular,h.irreg_pairs[0]);
                }
            }
        }
    }

    void Fill_Velocity_BC(const COMPONENT_LAYOUT_FEM<T>& cl)
    {
        typedef VECTOR<T,2> TV2;
        for(int b=0;b<cl.bc_v.m;b++)
        {
            const auto& bc=cl.bc_v(b);
            const auto& bl=cl.blocks(bc.b);
            MATRIX<T,2> M=bl.xform.M.Inverse();
            DOF_LAYOUT<TV> dl(cl,cl.reference_block_data(bl.ref_id),false);
            TV2 A=dl.cb->X(bc.bc_v.min_corner);
            TV2 B=dl.cb->X(bc.bc_v.max_corner-1);
            T width=(B-A).Magnitude();
            T k=bc.flow_rate*6/width;

            Visit_Dofs<false,true>(dl,LAYER_RANGE::ALL,bc.bc_v,bc.bc_e,
                [this,k,&bc,&M](const VISIT_ALL_DOFS<TV>& va)
                {
                    T z=(va.uv*((T)1-va.uv)).Product();
                    bc_v.Insert({bc.b,va.i},TV(M*-k*z*bc.normal));
                },
                [this,k,&bc,&M](const VISIT_ALL_DOFS<TV>& va)
                {
                    T z=(va.uv*((T)1-va.uv)).Product();
                    bc_e.Insert({bc.b,va.i},TV(M*-k*z*bc.normal));
                });
        }
    }

    void Build(const MATRIX_CONSTRUCTION_FEM<TV>& mc)
    {
        for(BLOCK_ID b(0);b<mc.cl.blocks.m;b++)
        {
            const auto& bl=mc.cl.blocks(b);
            const auto& rb=mc.cl.reference_block_data(bl.ref_id);
            DOF_LAYOUT<TV> dl(mc.cl,rb,false);
            first_v.Append(particles.Add_Elements(dl.counts.p));
            first_e.Append(element_edges.m*(TV::m-1)*3);
            Visit_Elements(dl,[this,&bl](const VISIT_ELEMENT_DATA<TV>& ve)
            {
                for(int i=0;i<ve.X.m;i++)
                    particles.X(first_v.Last()+ve.v(i))=xform(bl.xform,ve.X(i));
                mesh.elements.Append(ve.v+first_v.Last());
                element_edges.Append(ve.e+first_e.Last());
            });
            last_elem.Append(mesh.elements.m);
        }

        Fill_Global_Dof_Mapping(mc.cl);
        Fill_Velocity_BC(mc.cl);
    }

    TV Velocity(const TV& X) const
    {
        if(!tree)
            tree=new TREE(mesh,particles);

        return TV();
    }

    T Pressure(const TV& X) const
    {
        if(!tree)
            tree=new TREE(mesh,particles);

        return 0;
    }

    void Read(TYPED_ISTREAM input)
    {
    }

    void Write(TYPED_OSTREAM output) const
    {
    }
};
}
#endif
