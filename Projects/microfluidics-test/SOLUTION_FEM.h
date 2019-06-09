//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLUTION_FEM__
#define __SOLUTION_FEM__

#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
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

extern double comp_tol;

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

template<class T>
VECTOR<T,6> Velocity_Weights(const VECTOR<T,3>& x)
{
    T u=x(0),v=x(1);
    return {2*sqr(u)+4*u*v+2*sqr(v)-3*u-3*v+1,2*sqr(u)-u,2*sqr(v)-v,4*u*v,-4*u*v-4*sqr(v)+4*v,-4*sqr(u)-4*u*v+4*u};
}

template<class T>
VECTOR<T,3> Pressure_Weights(const VECTOR<T,3>& x)
{
    T u=x(0),v=x(1);
    return {1-u-v,u,v};
}

template<class T>
VECTOR<T,10> Velocity_Weights(const VECTOR<T,4>& x)
{
    T u=x(0),v=x(1),w=x(2);
    return {2*sqr(u)+4*u*v+4*u*w+2*sqr(v)+4*v*w+2*sqr(w)-3*u-3*v-3*w+1,
        2*sqr(u)-u,
        2*sqr(v)-v,
        2*sqr(w)-w,
        -4*sqr(u)-4*u*v-4*u*w+4*u,
        -4*u*v-4*sqr(v)-4*v*w+4*v,
        -4*u*w-4*v*w-4*sqr(w)+4*w,
        4*u*v,
        4*w*u,
        4*v*w};
}

template<class T>
VECTOR<T,4> Pressure_Weights(const VECTOR<T,4>& x)
{
    T u=x(0),v=x(1),w=x(2);
    return {1-u-v-w,u,v,w};
}

template<class TV>
struct SOLUTION_FEM
{
    typedef typename TV::SCALAR T;
    typedef typename MESH_POLICY<TV::m>::MESH MESH;
    typedef typename HIERARCHY_POLICY<TV>::TREE TREE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX SIMPLEX;

    GEOMETRY_PARTICLES<TV> particles;
    MESH mesh;
    ARRAY<VECTOR<int,(TV::m-1)*3> > element_edges;
    ARRAY<int,BLOCK_ID> first_v,first_e,last_elem;
    ARRAY<PAIR<BLOCK_ID,int> > dof_v,dof_e,dof_p;
    HASHTABLE<PAIR<BLOCK_ID,int>,TV> bc_v,bc_e;
    ARRAY<BLOCK_VECTOR<TV>,BLOCK_ID> sol_block_list;
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
        ARRAY<VECTOR<int,TV::m+1> > elements;
        for(BLOCK_ID b(0);b<mc.cl.blocks.m;b++)
        {
            const auto& bl=mc.cl.blocks(b);
            const auto& rb=mc.cl.reference_block_data(bl.ref_id);
            DOF_LAYOUT<TV> dl(mc.cl,rb,false);
            first_v.Append(particles.Add_Elements(dl.counts.p));
            first_e.Append(element_edges.m*(TV::m-1)*3);
            Visit_Elements(dl,[this,&bl,&elements](const VISIT_ELEMENT_DATA<TV>& ve)
            {
                for(int i=0;i<ve.X.m;i++)
                    particles.X(first_v.Last()+ve.v(i))=xform(bl.xform,ve.X(i));
                elements.Append(ve.v+first_v.Last());
                element_edges.Append(ve.e+first_e.Last());
            });
            last_elem.Append(elements.m);
        }
        mesh.Initialize_Mesh(particles.number,elements);

        Fill_Global_Dof_Mapping(mc.cl);
        Fill_Velocity_BC(mc.cl);
        sol_block_list=mc.rhs_block_list;
    }

    void Prepare_Hierarchy()
    {
        tree=new TREE(mesh,particles);
    }

    TV Velocity(const TV& X) const
    {
        auto getv=[this](int elem,int v)
        {
            const auto& dof=dof_v(v);
            if(dof.x<BLOCK_ID())
            {
                BLOCK_ID b=last_elem.Binary_Search(elem);
                const auto* bc=bc_v.Get_Pointer({b,v});
                if(bc) return *bc;
                else return TV();
            }
            else
                return sol_block_list(dof.x).Get_v(dof.y);
        };
        auto gete=[this](int elem,int e)
        {
            const auto& dof=dof_e(e);
            if(dof.x<BLOCK_ID())
            {
                BLOCK_ID b=last_elem.Binary_Search(elem);
                const auto* bc=bc_e.Get_Pointer({b,e});
                if(bc) return *bc;
                else return TV();
            }
            else
                return sol_block_list(dof.x).Get_e(dof.y);
        };

        ARRAY<int> hits;
        tree->Intersection_List(X,hits);
        for(int elem:hits)
        {
            auto simplex=SIMPLEX(particles.X.Subset(mesh.elements(elem)));
            VECTOR<T,TV::m+1> w=simplex.Barycentric_Coordinates(X);
            if(w.Min()>-comp_tol)
            {
                auto N=Velocity_Weights(w);
                VECTOR<TV,N.m> V;
                for(int i=0;i<mesh.elements(elem).m;i++)
                    V(i)=getv(elem,mesh.elements(elem)(i));
                for(int i=0;i<element_edges(elem).m;i++)
                    V(i+mesh.elements(elem).m)=gete(elem,element_edges(elem)(i));
                return V.Weighted_Sum(N);
            }
        }
        PHYSBAM_FATAL_ERROR();
    }

    T Pressure(const TV& X) const
    {
        auto get=[this](int elem,int v)
        {
            const auto& dof=dof_p(v);
            PHYSBAM_ASSERT(dof.x>=BLOCK_ID());
            return sol_block_list(dof.x).Get_p(dof.y);
        };
        ARRAY<int> hits;
        tree->Intersection_List(X,hits);
        for(int elem:hits)
        {
            auto simplex=SIMPLEX(particles.X.Subset(mesh.elements(elem)));
            VECTOR<T,TV::m+1> w=simplex.Barycentric_Coordinates(X);
            if(w.Min()>-comp_tol)
            {
                auto N=Pressure_Weights(w);
                VECTOR<T,N.m> P;
                for(int i=0;i<mesh.elements(elem).m;i++)
                    P(i)=get(elem,mesh.elements(elem)(i));
                return P.Weighted_Sum(N);
            }
        }
        PHYSBAM_FATAL_ERROR();
    }

    template<class RW> void Read(std::istream& input)
    {
        Read_Binary<RW>(input,particles,mesh,element_edges);
        Read_Binary<RW>(input,first_v,first_e,last_elem);
        Read_Binary<RW>(input,dof_v,dof_e,dof_p);
        Read_Binary<RW>(input,bc_v,bc_e);
        Read_Binary<RW>(input,sol_block_list);
    }

    template<class RW> void Write(std::ostream& output) const
    {
        Write_Binary<RW>(output,particles,mesh,element_edges);
        Write_Binary<RW>(output,first_v,first_e,last_elem);
        Write_Binary<RW>(output,dof_v,dof_e,dof_p);
        Write_Binary<RW>(output,bc_v,bc_e);
        Write_Binary<RW>(output,sol_block_list);
    }
};
}
#endif
