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
#include "ANALYTIC_FEM.h"
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
void Barycentric_Coordinates_Diff(const VECTOR<VECTOR<T,2>,3>& X,const VECTOR<T,2>& P,VECTOR<VECTOR<T,2>,3>& dw)
{
    typedef VECTOR<T,2> TV;
    TV u=X(1)-X(0),v=X(2)-X(0);
    T f=1/u.Cross(v)(0);
    dw(1)=TV(v(1),-v(0))*f;
    dw(2)=TV(-u(1),u(0))*f;
    dw(0)=-dw(1)-dw(2);
}

template<class T>
void Barycentric_Coordinates_Diff(const VECTOR<VECTOR<T,3>,4>& X,const VECTOR<T,3>& P,VECTOR<VECTOR<T,3>,4>& dw)
{
}

template<class T>
VECTOR<T,6> Velocity_Weights(const VECTOR<T,3>& x)
{
    T u=x(1),v=x(2);
    return {2*sqr(u)+4*u*v+2*sqr(v)-3*u-3*v+1,2*sqr(u)-u,2*sqr(v)-v,4*u*v,-4*u*v-4*sqr(v)+4*v,-4*sqr(u)-4*u*v+4*u};
}

template<class T>
MATRIX<T,2> Velocity_Gradient(const VECTOR<T,3>& w,const VECTOR<VECTOR<T,2>,3>& dw,const VECTOR<VECTOR<T,2>,6>& V)
{
    typedef VECTOR<T,2> TV;
    typedef MATRIX<T,2> TM;

    T u=w(1),v=w(2);
    VECTOR<TV,6> dNdw=
    {
        TV(4*u+4*v-3,4*u+4*v-3),
        TV(4*u-1,0),
        TV(0,4*v-1),
        TV(4*v,4*u),
        TV(-4*v,-4*u-8*v+4),
        TV(-8*u-4*v+4,-4*u)
    };
    TM dV;
    for(int i=0;i<6;i++)
    {
        VECTOR<TV,2> dwdx=dw.Remove_Index(0);
        TV dNdx=dwdx.Weighted_Sum(dNdw(i));
        dV+=TM::Outer_Product(V(i),dNdx);
    }
    return dV;
}

template<class T>
VECTOR<T,2> Pressure_Gradient(const VECTOR<T,3>& w,const VECTOR<VECTOR<T,2>,3>& dw,const VECTOR<T,3>& P)
{
    typedef VECTOR<T,2> TV;

    VECTOR<TV,3> dNdw=
    {
        TV(-1,-1),
        TV(1,0),
        TV(0,1),
    };
    TV dP;
    for(int i=0;i<3;i++)
    {
        VECTOR<TV,2> dwdx=dw.Remove_Index(0);
        TV dNdx=dwdx.Weighted_Sum(dNdw(i));
        dP+=P(i)*dNdx;
    }
    return dP;
}

template<class T>
VECTOR<T,3> Pressure_Gradient(const VECTOR<T,4>& w,const VECTOR<VECTOR<T,3>,4>& dw,const VECTOR<T,4>& P)
{
    typedef VECTOR<T,3> TV;
    return TV();
}

template<class T>
VECTOR<T,3> Pressure_Weights(const VECTOR<T,3>& x)
{
    T u=x(1),v=x(2);
    return {1-u-v,u,v};
}

template<class T>
VECTOR<T,10> Velocity_Weights(const VECTOR<T,4>& x)
{
    T u=x(1),v=x(2),w=x(3);
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
MATRIX<T,3> Velocity_Gradient(const VECTOR<T,4>& w,const VECTOR<VECTOR<T,3>,4>& dw,const VECTOR<VECTOR<T,3>,10>& V)
{
    typedef MATRIX<T,3> TM;
    return TM();
}

template<class T>
VECTOR<T,4> Pressure_Weights(const VECTOR<T,4>& x)
{
    T u=x(1),v=x(2),w=x(3);
    return {1-u-v-w,u,v,w};
}

template<class TV>
struct SOLUTION_FEM
{
    typedef typename TV::SCALAR T;
    typedef typename MESH_POLICY<TV::m>::MESH MESH;
    typedef typename HIERARCHY_POLICY<TV>::TREE TREE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX SIMPLEX;

    T intersection_tol=0;
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
            DOF_LAYOUT<TV> dl(cl,cl.reference_block_data(bl.ref_id),false);
            TV2 A=dl.cb->X(bc.bc_v.min_corner);
            TV2 B=dl.cb->X(bc.bc_v.max_corner-1);
            T width=(B-A).Magnitude();
            T k=bc.flow_rate*6/width;

            Visit_Dofs<false,true>(dl,LAYER_RANGE::ALL,bc.bc_v,bc.bc_e,
                [this,k,&bc](const VISIT_ALL_DOFS<TV>& va)
                {
                    T z=(va.uv*((T)1-va.uv)).Product();
                    bc_v.Insert({bc.b,first_v(bc.b)+va.i},TV(-k*z*bc.normal));
                },
                [this,k,&bc](const VISIT_ALL_DOFS<TV>& va)
                {
                    T z=(va.uv*((T)1-va.uv)).Product();
                    bc_e.Insert({bc.b,first_e(bc.b)+va.i},TV(-k*z*bc.normal));
                });
        }
    }

    void Fill_Analytic_Velocity_BC(const ANALYTIC_FEM<TV>& an)
    {
        T s=an.mc.cl.unit_s,m=an.mc.cl.unit_m;
        const COMPONENT_LAYOUT_FEM<T>& cl=an.mc.cl;
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            const auto& bl=cl.blocks(b);
            DOF_LAYOUT<TV> dl(cl,cl.reference_block_data(bl.ref_id),false);

            Visit_Dofs<true,false>(dl,LAYER_RANGE::ALL,bl.block->bc_v,bl.block->bc_e,
                [this,b,&bl,&an,s,m](const VISIT_ALL_DOFS<TV>& va)
                {
                    TV Z=xform(bl.xform,va.X);
                    bc_v.Insert({b,first_v(b)+va.i},TV(an.analytic_velocity->v(Z/m,0)*m/s));
                },
                [this,b,&bl,&an,s,m](const VISIT_ALL_DOFS<TV>& va)
                {
                    TV Z=xform(bl.xform,va.X);
                    bc_e.Insert({b,first_e(b)+va.i},TV(an.analytic_velocity->v(Z/m,0)*m/s));
                });
        }
    }

    void Build(const MATRIX_CONSTRUCTION_FEM<TV>& mc,const ANALYTIC_FEM<TV>* an)
    {
        intersection_tol=mc.cl.target_length/10;
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
        if(an)
            Fill_Analytic_Velocity_BC(*an);
        else
            Fill_Velocity_BC(mc.cl);
        sol_block_list=mc.rhs_block_list;
    }

    void Prepare_Hierarchy()
    {
        tree=new TREE(mesh,particles);
    }

    int Intersect(const TV& X) const
    {
        ARRAY<int> hits;
        tree->Intersection_List(X,hits,intersection_tol);
        for(int elem:hits)
        {
            auto simplex=SIMPLEX(particles.X.Subset(mesh.elements(elem)));
            VECTOR<T,TV::m+1> w=simplex.Barycentric_Coordinates(X);
            if(w.Min()>-comp_tol)
                return elem;
        }
        return -1;
    }

    TV Get_Vertex_Velocity(int elem,int v) const
    {
        const auto& dof=dof_v(v);
        if(dof.x<BLOCK_ID())
        {
            BLOCK_ID b=BLOCK_ID(std::upper_bound(last_elem.begin(),last_elem.end(),elem)-last_elem.begin());
            const auto* bc=bc_v.Get_Pointer({b,v});
            if(bc) return *bc;
            else return TV();
        }
        else
            return sol_block_list(dof.x).Get_v(dof.y);
    }

    TV Get_Edge_Velocity(int elem,int e) const
    {
        const auto& dof=dof_e(e);
        if(dof.x<BLOCK_ID())
        {
            BLOCK_ID b=BLOCK_ID(std::upper_bound(last_elem.begin(),last_elem.end(),elem)-last_elem.begin());
            const auto* bc=bc_e.Get_Pointer({b,e});
            if(bc) return *bc;
            else return TV();
        }
        else
            return sol_block_list(dof.x).Get_e(dof.y);
    }

    T Get_Pressure(int elem,int v) const
    {
        const auto& dof=dof_p(v);
        PHYSBAM_ASSERT(dof.x>=BLOCK_ID());
        return sol_block_list(dof.x).Get_p(dof.y);
    }

    TV Velocity(int elem,const TV& X,MATRIX<T,TV::m>* grad=0) const
    {
        auto simplex=SIMPLEX(particles.X.Subset(mesh.elements(elem)));
        VECTOR<T,TV::m+1> w=simplex.Barycentric_Coordinates(X);
        auto N=Velocity_Weights(w);
        VECTOR<TV,N.m> V;
        for(int i=0;i<mesh.elements(elem).m;i++)
            V(i)=Get_Vertex_Velocity(elem,mesh.elements(elem)(i));
        for(int i=0;i<element_edges(elem).m;i++)
            V(i+mesh.elements(elem).m)=Get_Edge_Velocity(elem,element_edges(elem)(i));
        if(grad)
        {
            VECTOR<TV,TV::m+1> dw;
            Barycentric_Coordinates_Diff(simplex.X,X,dw);
            *grad=Velocity_Gradient(w,dw,V);
        }
        return V.Weighted_Sum(N);
    }

    TV Velocity(const TV& X,MATRIX<T,TV::m>* grad=0) const
    {
        int elem=Intersect(X);
        PHYSBAM_ASSERT(elem>=0);
        return Velocity(elem,X,grad);
    }

    T Pressure(int elem,const TV& X,TV* grad=0) const
    {
        auto simplex=SIMPLEX(particles.X.Subset(mesh.elements(elem)));
        VECTOR<T,TV::m+1> w=simplex.Barycentric_Coordinates(X);
        auto N=Pressure_Weights(w);
        VECTOR<T,N.m> P;
        for(int i=0;i<mesh.elements(elem).m;i++)
            P(i)=Get_Pressure(elem,mesh.elements(elem)(i));
        if(grad)
        {
            VECTOR<TV,TV::m+1> dw;
            Barycentric_Coordinates_Diff(simplex.X,X,dw);
            *grad=Pressure_Gradient(w,dw,P);
        }
        return P.Weighted_Sum(N);
    }

    T Pressure(const TV& X,TV* grad=0) const
    {
        int elem=Intersect(X);
        PHYSBAM_ASSERT(elem>=0);
        return Pressure(elem,X,grad);
    }

    T Max_Pressure() const
    {
        T p=-FLT_MAX;
        for(const auto& bv:sol_block_list)
        {
            for(int i=0;i<bv.n.p;i++)
                p=std::max(p,bv.Get_p(i));
        }
        return p;
    }

    T Max_Velocity_Magnitude() const
    {
        T v=-FLT_MAX;
        for(const auto& bv:sol_block_list)
        {
            for(int i=0;i<bv.n.v;i++)
                v=std::max(v,bv.Get_v(i).Magnitude());
            for(int i=0;i<bv.n.e;i++)
                v=std::max(v,bv.Get_e(i).Magnitude());
        }
        for(const auto& bc:bc_v)
            v=std::max(v,bc.data.Magnitude());
        for(const auto& bc:bc_e)
            v=std::max(v,bc.data.Magnitude());
        return v;
    }

    template<class RW> void Read(std::istream& input)
    {
        Read_Binary<RW>(input,intersection_tol);
        Read_Binary<RW>(input,particles,mesh,element_edges);
        Read_Binary<RW>(input,first_v,first_e,last_elem);
        Read_Binary<RW>(input,dof_v,dof_e,dof_p);
        Read_Binary<RW>(input,bc_v,bc_e);
        Read_Binary<RW>(input,sol_block_list);
    }

    template<class RW> void Write(std::ostream& output) const
    {
        Write_Binary<RW>(output,intersection_tol);
        Write_Binary<RW>(output,particles,mesh,element_edges);
        Write_Binary<RW>(output,first_v,first_e,last_elem);
        Write_Binary<RW>(output,dof_v,dof_e,dof_p);
        Write_Binary<RW>(output,bc_v,bc_e);
        Write_Binary<RW>(output,sol_block_list);
    }
};
}
#endif
