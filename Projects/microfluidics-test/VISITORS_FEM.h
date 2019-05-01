//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VISITORS_FEM__
#define __VISITORS_FEM__
#include <Core/Math_Tools/INTERVAL.h>
#include "COMMON.h"
#include "COMPONENT_LAYOUT_FEM.h"

namespace PhysBAM{
// F(dof0,dof1,first_owns)
template<class F>
void Visit_Cross_Section_Dofs(INTERVAL<int> i0,INTERVAL<int> i1,bool uf0,bool uf1,bool m0,F func)
{
    int a=i1.min_corner,b=1,c=i0.min_corner,n=i0.Size();
    if(uf0==uf1)
    {
        a=i1.max_corner-1;
        b=-1;
    }
    if(uf0)
    {
        int m=n/2+(n&m0);
        for(int i=0;i<m;i++) func(c+i,a+b*i,true);
        for(int i=m;i<n;i++) func(c+i,a+b*i,false);
    }
    else
    {
        int m=n/2+(n&(1-m0));
        for(int i=0;i<m;i++) func(c+i,a+b*i,false);
        for(int i=m;i<n;i++) func(c+i,a+b*i,true);
    }
}

struct IRREGULAR_VISITOR
{
    BLOCK_ID b;
    int re,r0,r1;
    int ie,i0,i1;
    bool be,b0,b1; // regular owns e,v0,v1
    bool n0; // v0 not seen before for the irregular block
};

// F(const IRREGULAR_VISITOR& i)
template<class T,class F>
void Visit_Irregular_Cross_Section_Dofs(const ARRAY<BLOCK<T>,BLOCK_ID>& blocks,const IRREGULAR_CONNECTION& ic,F func)
{
    // edge: j=i*a+b
    // v0: j=i*a+c
    // v1: j=i*a+d
    auto& cs=blocks(ic.regular).block->cross_sections(ic.con_id);
    int a=1,b=0,c=0,d=1,n=ic.edge_on.m;
    if(!cs.own_first)
    {
        a=-1;
        b=n-1;
        d=n-1;
        c=n;
    }
    b+=cs.e.min_corner;
    c+=cs.v.min_corner;
    d+=cs.v.min_corner;

    int irreg_v0=n/2+1;
    int irreg_v1=n/2;
    BLOCK_ID last_b(-1);
    for(int i=0;i<n;i++)
    {
        bool b1=i<irreg_v1;
        bool b0=i<irreg_v0;
        auto& eo=ic.edge_on(i);
        IRREGULAR_VISITOR iv={eo.b,a*i+b,a*i+c,a*i+d,eo.e,eo.v0,eo.v1,b1,b0,b1,eo.b!=last_b};
        func(iv);
        last_b=eo.b;
    }
}

template<class CS,class FV,class FE>
void Visit_Regular_Cross_Section_Dofs(const CS& cs0,const CS& cs1,bool m0,FV func_v,FE func_e)
{
    Visit_Cross_Section_Dofs(cs0.v,cs1.v,cs0.own_first,cs1.own_first,m0,func_v);
    Visit_Cross_Section_Dofs(cs0.e,cs1.e,cs0.own_first,cs1.own_first,m0,func_e);
}

/*
RAW:

  2D counts: 
    v = [0,V) vertices
    e = [0,E) edges
    t = [0,T) triangles

  l = [0,L) layers
  k = [0,L] planes

  3D vertex dofs:
    v3 = v2 * (L+1) + k

  3D edge dofs:
    horizontal:
      e3 = e2 * (L+1) + k
    diagonal:
      e3 = e2 * L + E*(L+1) + l
    vertical:
      e3 = v2 * L + E*(2*L+1) + l

CONDENSED:

  3D vertex dofs (p):  [L+1]
    v3 = v2 * (L+1) + k

  3D vertex dofs (v):  [L-1]
    v3 = v2 * (L-1) + k

  3D edge dofs (e):
    horizontal:  [L-1]
      e3 = e2 * (L-1) + k
    diagonal:  [L]
      e3 = e2 * L + E*(L-1) + l
    vertical:  [L]
      e3 = v2 * L + E*(2*L-1) + l

 */


template<class TV> struct DOF_LAYOUT;

template<class T>
struct DOF_LAYOUT<VECTOR<T,2> >
{
    CANONICAL_BLOCK<T>* cb;
    DOF_COUNTS counts;

    DOF_LAYOUT(const COMPONENT_LAYOUT_FEM<T>& cl,const REFERENCE_BLOCK_DATA& rb,bool condensed)
        :cb(cl.blocks(rb.b).block)
    {
        counts=condensed?rb.num_dofs_d:rb.num_dofs_s;
    }
};

template<class T>
struct DOF_LAYOUT<VECTOR<T,3> >
{
    CANONICAL_BLOCK<T>* cb;
    DOF_COUNTS counts;
    int nl,nl1,d0,v0;
    T depth,dz;
    const ARRAY<int>& ticks_e;
    const ARRAY<int>& ticks_t;

    DOF_LAYOUT(const COMPONENT_LAYOUT_FEM<T>& cl,const REFERENCE_BLOCK_DATA& rb,bool condensed)
        :cb(cl.blocks(rb.b).block),nl(cl.depth_layers),depth(cl.depth),dz(depth/nl),
        ticks_e(rb.ticks_e),ticks_t(rb.ticks_t)
    {
        nl1=nl+1-2*condensed;
        counts=condensed?rb.num_dofs_d:rb.num_dofs_s;
        int num_eh=counts.e*nl1,num_ed=counts.e*nl,num_ev=counts.v*nl;
        d0=num_eh;
        v0=num_eh+num_ed;
        counts.e=num_eh+num_ed+num_ev;
        counts.v*=nl1;
        counts.p*=nl+1;
    }

    int Vertex(int v, int l) const {return v * nl1 + l;}
    int Vertex_p(int v, int l) const {return v * (nl+1) + l;}
    int Edge_h(int e, int l) const {return e * nl1 + l;}
    int Edge_d(int e, int l) const {return e * nl + d0 + l;}
    int Edge_v(int v, int l) const {return v * nl + v0 + l;}

    VECTOR<T,3> X(int v, int l) const {return cb->X(v).Append(l*dz);}
};


static const int move_order_table[7][3]=
{
    {-7,-7,-7},
    {1,0,2},
    {2,1,0},
    {1,2,0},
    {0,2,1},
    {0,1,2},
    {2,0,1}
};

template<class TV>
struct VISIT_ELEMENT_DATA
{
    int tri; // template triangle (2D)
    VECTOR<int,TV::m+1> v; // vertex dofs (3D)
    VECTOR<int,3*(TV::m-1)> e; // edge dofs (3D)
    VECTOR<TV,TV::m+1> X;
};

// func(const VISIT_ELEMENT_DATA<T>& vt);
template<class T,class F>
void Visit_Elements(const DOF_LAYOUT<VECTOR<T,2> >& dl,F func)
{
    typedef VECTOR<T,2> TV;
    int m=dl.cb->E.m;
    for(int tri=0;tri<m;tri++)
    {
        VISIT_ELEMENT_DATA<TV> vt;
        vt.tri=tri;
        vt.v=dl.cb->E(tri);
        for(int i=0;i<3;i++)
            vt.e(i)=dl.cb->element_edges(tri)(i).x;
        vt.X={dl.cb->X(vt.v(0)),dl.cb->X(vt.v(1)),dl.cb->X(vt.v(2))};
        func(vt);
    }
}
        
// func(const VISIT_ELEMENT_DATA<T>& vt);
template<class T,class F>
void Visit_Elements(const DOF_LAYOUT<VECTOR<T,3> >& dl,F func)
{
    typedef VECTOR<T,3> TV;
    int m=dl.cb->E.m;

    for(int tri=0;tri<m;tri++)
    {
        VISIT_ELEMENT_DATA<TV> vt[3],mirror_vt[3];
        VECTOR<int,3> P=dl.cb->E(tri),S;
        for(int i=0;i<3;i++) S(i)=dl.cb->element_edges(tri)(i).x;
    
        auto fill_vt=[=](auto& keys,VISIT_ELEMENT_DATA<TV>& d)
            {
                d.tri=tri;
                for(int i=0;i<4;i++)
                {
                    d.v(i)=dl.Vertex(P(keys(i).x),keys(i).y);
                    d.X(i)=dl.X(P(keys(i).x),keys(i).y);
                }
                for(int i=0,e=0;i<4;i++)
                    for(int j=i+1;j<4;j++)
                    {
                        int r=0;
                        if(keys(i).x==keys(j).x) r=dl.Edge_v(P(keys(i).x),0);
                        else
                        {
                            int s=S(3-keys(i).x-keys(j).x);
                            if(keys(i).y==keys(j).y) r=dl.Edge_h(s,keys(i).y);
                            else r=dl.Edge_d(s,0);
                        }
                        d.e(e++)=r;
                    }
            };

        VECTOR<VECTOR<int,2>,4> keys({0,0},{1,0},{2,0},{0,1});
        const int* move_order=move_order_table[dl.ticks_t(tri)];
        for(int i=0;i<3;i++)
        {
            keys(3).x=move_order[i];
            fill_vt(keys,vt[i]);
            keys(move_order[i])=keys(3);
        }
        keys=VECTOR<VECTOR<int,2>,4>({0,0},{1,0},{2,0},{0,1});
        for(int i=0;i<3;i++)
        {
            int j=2-i;
            keys(3).x=move_order[j];
            fill_vt(keys,mirror_vt[i]);
            keys(move_order[j])=keys(3);
            mirror_vt[i].v+=dl.nl-1;
            mirror_vt[i].e+=dl.nl-1;
            mirror_vt[i].X+=TV(0,0,(dl.nl-1)*dl.dz);
        }

        for(int l=0;l<dl.nl;l++)
        {
            VISIT_ELEMENT_DATA<TV>* vtp=vt;
            if(l==dl.nl-1) vtp=mirror_vt;
            for(int i=0;i<3;i++)
            {
                func(vtp[i]);
                vtp[i].v+=1;
                vtp[i].e+=1;
                vtp[i].X+=TV(0,0,dl.dz);
            }
        }
    }
}

template<class TV>
struct VISIT_FACE_DATA
{
    int edge; // template edge (2D)
    VECTOR<int,TV::m> v; // vertex dofs (3D)
    VECTOR<int,2*TV::m-3> e; // edge dofs (3D)
    VECTOR<TV,TV::m> X;
};

// func(const VISIT_FACE_DATA<T>& vt)
template<class T,class F>
void Visit_Faces(const DOF_LAYOUT<VECTOR<T,2> >& dl,INTERVAL<int> bc_e,F func)
{
    typedef VECTOR<T,2> TV;

    for(int edge:bc_e)
    {
        VISIT_FACE_DATA<TV> vt;
        vt.edge=edge;
        vt.v=dl.cb->S(edge);
        vt.e(0)=edge;
        vt.X={dl.cb->X(vt.v(0)),dl.cb->X(vt.v(1))};
        func(vt);
    }
}

// func(const VISIT_FACE_DATA<T>& vt)
template<class T,class F>
void Visit_Faces(const DOF_LAYOUT<VECTOR<T,3> >& dl,INTERVAL<int> bc_e,F func)
{
    typedef VECTOR<T,3> TV;

    auto fill_vt=[=](int edge,VISIT_FACE_DATA<TV>* vt,bool mirrored)
    {
        VECTOR<int,2> P=dl.cb->S(edge);
        if((dl.ticks_e(edge) && !mirrored) || (!dl.ticks_e(edge) && mirrored)) std::swap(P.x,P.y);
        vt[0].edge=edge;
        vt[0].v={dl.Vertex(P.x,0),dl.Vertex(P.y,0),dl.Vertex(P.y,1)};
        vt[0].e(0)=dl.Edge_v(P.y,0);
        vt[0].e(1)=dl.Edge_d(edge,0);
        vt[0].e(2)=dl.Edge_h(edge,0);
        vt[0].X={dl.X(P.x,0),dl.X(P.y,0),dl.X(P.y,1)};

        vt[1].edge=edge;
        vt[1].v={dl.Vertex(P.x,1),dl.Vertex(P.x,0),dl.Vertex(P.y,1)};
        vt[1].e(0)=dl.Edge_d(edge,0);
        vt[1].e(1)=dl.Edge_h(edge,1);
        vt[1].e(2)=dl.Edge_v(P.x,0);
        vt[1].X={dl.X(P.x,1),dl.X(P.x,0),dl.X(P.y,1)};
    };

    for(int edge:bc_e)
    {
        VISIT_FACE_DATA<TV> vt[2],mirror_vt[2];
        fill_vt(edge,vt,false);
        fill_vt(edge,mirror_vt,true);
        for(int i=0;i<2;i++)
        {
            mirror_vt[i].v+=dl.nl-1;
            mirror_vt[i].e+=dl.nl-1;
            mirror_vt[i].X+=TV(0,0,(dl.nl-1)*dl.dz);
        }

        for(int l=0;l<dl.nl;l++)
        {
            VISIT_FACE_DATA<TV>* vtp=vt;
            if(l==dl.nl-1) vtp=mirror_vt;
            for(int i=0;i<2;i++)
            {
                func(vt[i]);
                vtp[i].v+=1;
                vtp[i].e+=1;
                vtp[i].X+=TV(0,0,dl.dz);
            }
        }
    }
}

template<class TV>
struct VISIT_ALL_DOFS
{
    typedef typename TV::SCALAR T;
    int i;
    TV X;
    VECTOR<T,TV::m-1> uv;
};

// fv(int v,const TV& X,VECTOR<T,1>(u)); u=[0,1]
// fe(int e,const TV& X,VECTOR<T,1>(u)); u=[0,1]
template<bool use_X,bool use_uv,class T,class AR,class FV,class FE>
void Visit_Dofs(const DOF_LAYOUT<VECTOR<T,2> >& dl,const AR& bc_v,const AR& bc_e,FV fv,FE fe)
{
    typedef VECTOR<T,2> TV;

    TV A,u;
    if(use_uv)
    {
        A=dl.cb->X(*bc_v.begin());
        auto it=bc_v.end();
        u=dl.cb->X(*--it)-A;
        u/=u.Magnitude_Squared();
    }
    
    for(int v:bc_v)
    {
        VISIT_ALL_DOFS<TV> va;
        va.i=v;
        TV X;
        if(use_X || use_uv) X=dl.cb->X(v);
        if(use_X) va.X=X;
        if(use_uv) va.uv.x=(X-A).Dot(u);
        fv(va);
    }
    for(int e:bc_e)
    {
        VISIT_ALL_DOFS<TV> va;
        va.i=e;
        TV X=dl.cb->X.Subset(dl.cb->S(e)).Sum()/2;
        if(use_X) va.X=X;
        if(use_uv) va.uv.x=(X-A).Dot(u);
        fe(va);
    }
}

// fv(int v,const TV& X,VECTOR<T,2>(u,z)); u=[0,1], z=[0,1]
// fe(int e,const TV& X,VECTOR<T,2>(u,z)); u=[0,1], z=[0,1]
template<bool use_X,bool use_uv,class T,class AR,class FV,class FE>
void Visit_Dofs(const DOF_LAYOUT<VECTOR<T,3> >& dl,const AR& bc_v,const AR& bc_e,FV fv,FE fe)
{
    typedef VECTOR<T,2> TV2;
    typedef VECTOR<T,3> TV;

    TV2 A,u;
    if(use_uv)
    {
        A=dl.cb->X(*bc_v.begin());
        auto it=bc_v.end();
        u=dl.cb->X(*--it)-A;
        u/=u.Magnitude_Squared();
    }
    
    for(int v:bc_v)
    {
        VISIT_ALL_DOFS<TV> va;
        TV2 X=dl.cb->X(v);
        if(use_X) va.X=X.Append(0);
        if(use_uv) va.uv.x=(X-A).Dot(u);
        int v0=dl.Vertex(v,0);
        for(int i=0;i<dl.nl+1;i++)
        {
            va.i=v0+i;
            if(use_uv) va.uv.y=(T)i/dl.nl;
            if(use_X) va.X.z=i*dl.dz;
            fv(va);
        }
        int e0=dl.Edge_v(v,0);
        for(int i=0;i<dl.nl;i++)
        {
            va.i=e0+i;
            if(use_uv) va.uv.y=(i+(T).5)/dl.nl;
            if(use_X) va.X.z=(i+(T).5)*dl.dz;
            fe(va);
        }
    }
    for(int e:bc_e)
    {
        VISIT_ALL_DOFS<TV> va;
        TV2 X=dl.cb->X.Subset(dl.cb->S(e)).Sum()/2;
        if(use_X) va.X=X.Append(0);
        if(use_uv) va.uv.x=(X-A).Dot(u);
        int e0=dl.Edge_h(e,0);
        for(int i=0;i<dl.nl+1;i++)
        {
            va.i=e0+i;
            if(use_uv) va.uv.y=(T)i/dl.nl;
            if(use_X) va.X.z=i*dl.dz;
            fe(va);
        }
        int e1=dl.Edge_d(e,0);
        for(int i=0;i<dl.nl;i++)
        {
            va.i=e1+i;
            if(use_uv) va.uv.y=(i+(T).5)/dl.nl;
            if(use_X) va.X.z=(i+(T).5)*dl.dz;
            fe(va);
        }
    }
}

// func(int dest,int src)
template<class T,class FV,class FE,class FP>
void Visit_Dof_Pairs(const DOF_LAYOUT<VECTOR<T,2> >& dl0,const DOF_LAYOUT<VECTOR<T,2> >& dl1,const DOF_PAIRS& dp,FV fv,FE fe,FP fp)
{
    for(auto p:dp.v) fv(p.x,p.y);
    for(auto p:dp.e) fe(p.x,p.y);
    for(auto p:dp.p) fp(p.x,p.y);
}

// func(int dest,int src)
template<class T,class FV,class FE,class FP>
void Visit_Dof_Pairs(const DOF_LAYOUT<VECTOR<T,3> >& dl0,const DOF_LAYOUT<VECTOR<T,3> >& dl1,const DOF_PAIRS& dp,FV fv,FE fe,FP fp)
{
    int nl=dl0.nl;
    int nl1=dl0.nl1;
    assert(nl==dl1.nl);
    assert(nl1==dl1.nl-1);
    for(auto p:dp.v)
    {
        int v0=dl0.Vertex(p.x,0);
        int v1=dl1.Vertex(p.y,0);
        for(int i=0;i<nl1;i++) fv(v0+i,v1+i+1);
        int e0=dl0.Edge_v(p.x,0);
        int e1=dl1.Edge_v(p.y,0);
        for(int i=0;i<nl;i++) fe(e0+i,e1+i);
    }
    for(auto p:dp.e)
    {
        int h0=dl0.Edge_h(p.x,0);
        int h1=dl1.Edge_h(p.y,0);
        for(int i=0;i<nl1;i++) fe(h0+i,h1+i+1);
        int d0=dl0.Edge_d(p.x,0);
        int d1=dl1.Edge_d(p.y,0);
        for(int i=0;i<nl;i++) fe(d0+i,d1+i);
    }
    for(auto p:dp.p)
    {
        int v0=dl0.Vertex_p(p.x,0);
        int v1=dl1.Vertex_p(p.y,0);
        for(int i=0;i<nl+1;i++) fp(v0+i,v1+i);
    }
}

// func(int i,const TV& X)
template<class TV,class FV,class FE>
void Visit_Dofs(const DOF_LAYOUT<TV>& dl,FV fv,FE fe)
{
    Visit_Dofs<true,false>(dl,INTERVAL<int>(0,dl.cb->X.m),INTERVAL<int>(0,dl.cb->S.m),fv,fe);
}

template<class T,class FV,class FE,class FP>
void Visit_Compressed_Dofs(const DOF_LAYOUT<VECTOR<T,2> >& dl,const REFERENCE_BLOCK_DATA& rd,FV fv,FE fe,FP fp)
{
    typedef VECTOR<T,2> TV;

    for(int i=0;i<dl.cb->X.m;i++)
    {
        int j=rd.dof_map_v(i);
        if(j>=0) fv(j,dl.cb->X(i));
    }

    for(int i=0;i<dl.cb->S.m;i++)
    {
        int j=rd.dof_map_e(i);
        if(j>=0) fe(j,dl.cb->X.Subset(dl.cb->S(i)).Sum()/2);
    }

    for(int i=0;i<dl.cb->X.m;i++)
    {
        int j=rd.dof_map_p(i);
        if(j>=0) fp(j,dl.cb->X(i));
    }
}

template<class T,class FV,class FE,class FP>
void Visit_Compressed_Dofs(const DOF_LAYOUT<VECTOR<T,3> >& dl,const REFERENCE_BLOCK_DATA& rd,FV fv,FE fe,FP fp)
{
    typedef VECTOR<T,2> TV2;

    for(int i=0;i<dl.cb->X.m;i++)
    {
        int j=rd.dof_map_v(i);
        if(j>=0)
        {
            TV2 X=dl.cb->X(i);
            int v0=dl.Vertex(j,0);
            for(int k=0;k<dl.nl1;k++) fv(v0+k,X.Append((T)(k+1)*dl.dz));
            int e0=dl.Edge_v(j,0);
            for(int k=0;k<dl.nl;k++) fe(e0+k,X.Append((k+(T).5)*dl.dz));
        }
    }

    for(int i=0;i<dl.cb->S.m;i++)
    {
        int j=rd.dof_map_e(i);
        if(j>=0)
        {
            TV2 X=dl.cb->X.Subset(dl.cb->S(i)).Sum()/2;
            int e0=dl.Edge_h(j,0);
            for(int k=0;k<dl.nl1;k++) fe(e0+k,X.Append((T)(k+1)*dl.dz));
            int e1=dl.Edge_d(j,0);
            for(int k=0;k<dl.nl;k++) fe(e1+k,X.Append((k+(T).5)*dl.dz));
        }
    }

    for(int i=0;i<dl.cb->X.m;i++)
    {
        int j=rd.dof_map_p(i);
        if(j>=0)
        {
            TV2 X=dl.cb->X(i);
            int v0=dl.Vertex_p(j,0);
            for(int k=0;k<dl.nl+1;k++) fp(v0+k,X.Append((T)k*dl.dz));
        }
    }
}

}
#endif
