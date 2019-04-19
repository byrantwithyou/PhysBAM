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


 */

static const int move_order_table[7][3]={};

struct VISIT_TET_DATA
{
    int tri; // template triangle (2D)
    VECTOR<int,4> v; // vertex dofs (3D)
    VECTOR<int,6> e; // edge dofs (3D)
};

// func(const VISIT_TET_DATA& vt);
template<class T,class F>
void Visit_Triangle_Tetrahedra(const COMPONENT_LAYOUT_FEM<T>& cl,const REFERENCE_BLOCK_DATA& rb,int tri,F func)
{
    int* move_order=move_order_table[rb.ticks_t(tri)];
    auto* cb=cl.blocks(rb.b).block;
    int nl=cl.depth_layers,ne=cb->S.m,np=cb->X.m;
    
    VISIT_TET_DATA vt[3];
    VECTOR<int,3> P=cb->E(tri),S;
    for(int i=0;i<3;i++) S(i)=cb->element_edges(i).x;
    
    auto fill_vt=[=](auto& keys,VISIT_TET_DATA& d)
        {
            d.tri=tri;
            for(int i=0;i<4;i++)
                d.v(i)=P(keys(i).x)*(nl+1)+keys(i).y;
            for(int i=0,e=0;i<4;i++)
                for(int j=i+1;j<4;j++)
                {
                    int r=0;
                    if(keys(i).x==keys(j).x) r=P(keys(i).x)*nl+ne*(2*nl+1);
                    else
                    {
                        int s=S(3-keys(i).x-keys(j).x);
                        if(keys(i).y==keys(j).y) r=s*(nl+1);
                        else r=s*nl+ne*(nl+1);
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

    for(int l=0;l<nl;l++)
        for(int i=0;i<3;i++)
        {
            func(vt[i]);
            vt[i].v+=1;
            vt[i].e+=1;
        }
}

template<class T,class F>
void Visit_Block_Tetrahedra(const COMPONENT_LAYOUT_FEM<T>& cl,const REFERENCE_BLOCK_DATA& rb,F func)
{
    int m=cl.blocks(rb.b).block->E.m;
    for(int i=0;i<m;i++)
        Visit_Triangle_Tetrahedra(cl,rb,i,func);
}

template<class T,class F>
void Visit_Edge_Faces(const COMPONENT_LAYOUT_FEM<T>& cl,const REFERENCE_BLOCK_DATA& rb,int edge,F func)
{
    
}

template<class T,class F>
void Visit_Block_Faces(const COMPONENT_LAYOUT_FEM<T>& cl,const REFERENCE_BLOCK_DATA& rb,F func)
{
    int m=cl.blocks(rb.b).block->S.m;
    for(int i=0;i<m;i++)
        Visit_Edge_Faces(cl,rb,i,func);
}


}
#endif
