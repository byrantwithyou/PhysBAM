//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VISITORS_FEM__
#define __VISITORS_FEM__

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

    DOF_LAYOUT(const COMPONENT_LAYOUT_FEM<T>& cl,const REFERENCE_BLOCK_DATA& rb)
        :cb(cl.blocks(rb.b).block)
    {}
};

template<class T>
struct DOF_LAYOUT<VECTOR<T,3> >
{
    CANONICAL_BLOCK<T>* cb;
    int nl,nl1,d0,v0;
    int num_v,num_e,num_p;
    T depth,dz;
    const ARRAY<int>& ticks_e;
    const ARRAY<int>& ticks_t;

    DOF_LAYOUT(const COMPONENT_LAYOUT_FEM<T>& cl,const REFERENCE_BLOCK_DATA& rb,int num_edges,bool condensed)
        :cb(cl.blocks(rb.b).block),nl(cl.depth_layers),depth(cl.depth),dz(depth/nl),
        ticks_e(rb.ticks_e),ticks_t(rb.ticks_t)
    {
        nl1=nl+1-2*condensed;
        num_v=nl1;
        num_e=condensed?rb.num_dofs_d.e:rb.num_dofs_s.e;
        num_p=nl+1;
        d0=num_e*nl1;
        v0=d0+num_e*nl;
    }

    int Vertex(int v, int l, bool is_p=false) const {return v * (nl1+2*is_p) + l;}
    int Edge_h(int e, int l) const {return e * nl1 + l;}
    int Edge_d(int e, int l) const {return e * nl + d0 + l;}
    int Edge_v(int v, int l) const {return v * nl + v0 + l;}

    VECTOR<T,3> X(int v, int l) const {return cb->X(v).Append(l*dz);}
};


static const int move_order_table[7][3]={}; // TODO

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
void Visit_Elements(const DOF_LAYOUT<VECTOR<T,3> >& dl,F func)
{
    typedef VECTOR<T,3> TV;
    int m=dl.cb->E.m;
    for(int tri=0;tri<m;tri++)
    {
        int* move_order=move_order_table[dl.ticks_t(tri)];
    
        VISIT_ELEMENT_DATA<T> vt[3];
        VECTOR<int,3> P=dl.cb->E(tri),S;
        for(int i=0;i<3;i++) S(i)=dl.cb->element_edges(i).x;
    
        auto fill_vt=[=](auto& keys,VISIT_ELEMENT_DATA<T>& d)
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
        for(int i=0;i<3;i++)
        {
            keys(3).x=move_order[i];
            fill_vt(keys,vt[i]);
            keys(move_order[i])=keys(3);
        }

        for(int l=0;l<dl.nl;l++)
            for(int i=0;i<3;i++)
            {
                func(vt[i]);
                vt[i].v+=1;
                vt[i].e+=1;
                vt[i].X+=TV(0,0,dl.dz);
            }
    }
}

template<class TV>
struct VISIT_FACE_DATA
{
    int edge; // template edge (2D)
    VECTOR<int,TV::m> v; // vertex dofs (3D)
    VECTOR<int,2*TV::m-1> e; // edge dofs (3D)
    VECTOR<TV,TV::m> X;
};

// func(const VISIT_FACE_DATA<T>& vt)
template<class T,class F>
void Visit_Faces(const DOF_LAYOUT<VECTOR<T,3> >& dl,F func)
{
    typedef VECTOR<T,3> TV;

    int m=dl.cb->S.m;
    for(int edge=0;edge<m;edge++)
    {
        VECTOR<int,2> P=dl.cb->S(edge);
        if(dl.ticks_e(edge)) std::swap(P.x,P.y);
        VISIT_FACE_DATA<T> vt[2];
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

        for(int l=0;l<dl.nl;l++)
            for(int i=0;i<2;i++)
            {
                func(vt[i]);
                vt[i].v+=1;
                vt[i].e+=1;
                vt[i].X+=TV(0,0,dl.dz);
            }
    }
}

// func(int dest,int src)
template<class T,class FV,class FE,class FP>
void Visit_Dof_Pairs(const DOF_LAYOUT<VECTOR<T,2> >& dl,const DOF_PAIRS& dp,FV fv,FE fe,FP fp)
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
    assert(nl1==dl1.nl1);
    for(auto p:dp.v)
    {
        int v0=dl0.Vertex(p.x,0,false);
        int v1=dl1.Vertex(p.y,0,false);
        for(int i=0;i<nl1;i++) fv(v0+i,v1+i);
        int e0=dl0.Edge_v(p.x,0);
        int e1=dl1.Edge_v(p.y,0);
        for(int i=0;i<nl;i++) fe(e0+i,e1+i);
    }
    for(auto p:dp.e)
    {
        int h0=dl0.Edge_h(p.x,0);
        int h1=dl1.Edge_h(p.y,0);
        for(int i=0;i<nl1;i++) fe(h0+i,h1+i);
        int d0=dl0.Edge_d(p.x,0);
        int d1=dl1.Edge_d(p.y,0);
        for(int i=0;i<nl;i++) fe(d0+i,d1+i);
    }
    for(auto p:dp.p)
    {
        int v0=dl0.Vertex(p.x,0,true);
        int v1=dl1.Vertex(p.y,0,true);
        for(int i=0;i<nl+1;i++) fp(v0+i,v1+i);
    }
}

}
#endif
