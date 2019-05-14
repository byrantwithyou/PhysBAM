//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
//#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Data_Structures/TUPLE.h>
#include <Core/Data_Structures/UNION_FIND.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Math_Tools/RANGE.h>
//#include <Core/Matrices/MATRIX_MXN.h>
//#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
//#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
//#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
//#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
//#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <fstream>
//#include "CACHED_ELIMINATION_MATRIX.h"
#include "COMPONENT_BC.h"
#include "COMPONENT_CHANGE.h"
#include "COMPONENT_JOINT.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include "COMPONENT_PIPE.h"
#include "VISITORS_FEM.h"
#include <tuple>
namespace PhysBAM{

double comp_tol=1e-10;

//#####################################################################
// Function Destructor
//#####################################################################
template<class T> COMPONENT_LAYOUT_FEM<T>::
~COMPONENT_LAYOUT_FEM()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Update_Masters()
{
    for(BLOCK_ID b0(0);b0<blocks.m;b0++)
    {
        for(CON_ID c0(0);c0<blocks(b0).connections.m;c0++)
        {
            auto& n0=blocks(b0).connections(c0);
            BLOCK_ID b1=n0.id;
            CON_ID c1=n0.con_id;
            if(!n0.is_regular)
            {
                n0.master=true;
                continue;
            }
            if(b0>b1) continue; // smaller index does the logic
            auto& n1=blocks(b1).connections(c1);
            bool m=true;
            if(c0==CON_ID())
            {
                if(c1!=CON_ID()) m=true;
                else
                {
                    if(blocks(b0).connections.m==CON_ID(2)) m=true;
                    else if(blocks(b1).connections.m==CON_ID(2)) m=false;
                }
            }
            else
            {
                if(!c1) m=false;
                else if(blocks(b0).connections.m==CON_ID(2)) m=true;
                else if(blocks(b1).connections.m==CON_ID(2)) m=false;
            }
            n0.master=m;
            n1.master=!m;
        }
    }
}
//#####################################################################
// Function Separates_Dofs
//#####################################################################
template<class T> int COMPONENT_LAYOUT_FEM<T>::
Separates_Dofs(BLOCK_ID b)
{
    CANONICAL_BLOCK<T>* cb=blocks(b).block;
    int mask=0;
    for(CON_ID c(0);c<blocks(b).connections.m;c++)
        if(blocks(b).connections(c).master)
            mask|=1<<Value(c);
    auto pr=separates_dofs.Insert({cb,mask},0);
    if(!pr.y) return *pr.x;

    int valid=0;
    ARRAY<CON_ID> owner(cb->X.m,use_init,CON_ID(-1));
    for(CON_ID i(0);i<cb->cross_sections.m;i++)
    {
        const auto& cs=cb->cross_sections(i);
        int master=(mask>>Value(i))&1;
        INTERVAL<int> iv=cs.v;
        int m=(iv.min_corner+iv.max_corner)/2;
        if(cs.own_first) iv.min_corner=m+(iv.Size()&master);
        else iv.max_corner=m+(iv.Size()&(1-master));
        for(int j:iv) owner(j)=i;
    }
    for(auto e:cb->E)
    {
        VECTOR<CON_ID,3> f(owner.Subset(e));
        CON_ID o(-1);
        for(auto i:f)
        {
            if(i==CON_ID(-1) || i==o) continue;
            if(o==CON_ID(-1))
            {
                o=i;
                continue;
            }
            valid|=1<<Value(o);
            valid|=1<<Value(i);
        }
    }
    *pr.x=valid;
    return valid;
}
int Merge_Dofs(ARRAY<int>& index_map,int num_dofs0,int num_dofs1,
    INTERVAL<int> i0,INTERVAL<int> i1,bool uf0,bool uf1)
{
    index_map.Resize(num_dofs1,init_all,-1);

    Visit_Cross_Section_Dofs(i0,i1,uf0,uf1,true,[&index_map](int a,int b,bool){index_map(b)=a;});

    int k=num_dofs0;
    for(int i=0;i<index_map.m;i++)
        if(index_map(i)<0)
            index_map(i)=k++;
    return k;
}
PAIR<INTERVAL<int>,bool> Map_Interval(INTERVAL<int> iv,const ARRAY<int>& index_map)
{
    int a=index_map(iv.min_corner);
    int b=index_map(iv.max_corner-1);
    if(a<=b) return {{a,b+1},false};
    return {{b,a+1},true};
}
template<class CS>
CS Map_Cross_Section(CS cs,const ARRAY<int>& index_v_map,const ARRAY<int>& index_e_map)
{
    auto pv=Map_Interval(cs.v,index_v_map);
    auto pe=Map_Interval(cs.e,index_e_map);
    assert(pv.y==pe.y);
    cs.v=pv.x;
    cs.e=pe.x;
    if(pv.y) cs.own_first=!cs.own_first;
    return cs;
}
//#####################################################################
// Function Merge_Canonical_Blocks
//#####################################################################
template<class T> TRIPLE<CANONICAL_BLOCK<T>*,ARRAY<int>,ARRAY<int> >& COMPONENT_LAYOUT_FEM<T>::
Merge_Canonical_Blocks(CANONICAL_BLOCK<T>* cb0,CON_ID con_id0,XFORM<TV> xf0,
    CANONICAL_BLOCK<T>* cb1,CON_ID con_id1,XFORM<TV> xf1)
{
    auto pr=merge_canonical_blocks.Insert(std::make_tuple(cb0,con_id0,cb1,con_id1),{});
    if(!pr.y) return *pr.x;

    const CROSS_SECTION& cs0=cb0->cross_sections(con_id0);
    const CROSS_SECTION& cs1=cb1->cross_sections(con_id1);

    ARRAY<int>& index_v_map=pr.x->y;
    ARRAY<int>& index_e_map=pr.x->z;
    int num_v_dofs=Merge_Dofs(index_v_map,cb0->X.m,cb1->X.m,cs0.v,cs1.v,
        cs0.own_first,cs1.own_first);
    int num_e_dofs=Merge_Dofs(index_e_map,cb0->S.m,cb1->S.m,cs0.e,cs1.e,
        cs0.own_first,cs1.own_first);

    pr.x->x=new CANONICAL_BLOCK<T>;
    CANONICAL_BLOCK<T>* cb=pr.x->x;
    for(CON_ID i(0);i<cb0->cross_sections.m;i++)
        if(i!=con_id0)
            cb->cross_sections.Append(cb0->cross_sections(i));
    for(CON_ID i(0);i<cb1->cross_sections.m;i++)
        if(i!=con_id1)
            cb->cross_sections.Append(
                Map_Cross_Section(cb1->cross_sections(i),index_v_map,index_e_map));

    XFORM<TV> M01i=xf0.Inverse()*xf1;

    cb->X=cb0->X;
    cb->X.Resize(num_v_dofs);
    for(int i=0;i<cb1->X.m;i++)
    {
        TV X=M01i*cb1->X(i);
        int j=index_v_map(i);
        if(j>=cb0->X.m) cb->X(j)=X;
        else assert((cb->X(j)-X).Magnitude()<(T)1e-6);
    }

    cb->S=cb0->S;
    cb->S.Resize(num_e_dofs);
    cb->ticks=cb0->ticks;
    cb->ticks.Resize(num_e_dofs);
    for(int i=0;i<cb1->S.m;i++)
    {
        int j=index_e_map(i);
        IV s(index_v_map.Subset(cb1->S(i)));
        if(j>=cb0->S.m)
        {
            cb->S(j)=s;
            cb->ticks(j)=cb1->ticks(i);
        }
        else
        {
            auto s0=cb->S(j);
            s.Sort();
            s0.Sort();
            assert(s0==s);
        }
    }

    cb->bc_v=cb0->bc_v;
    cb->bc_e=cb0->bc_e;
    for(int v:cb1->bc_v)
    {
        int j=index_v_map(v);
        if(j>=cb0->X.m) cb->bc_v.Append(j);
    }
    for(int e:cb1->bc_e)
    {
        int j=index_e_map(e);
        if(j>=cb0->S.m) cb->bc_e.Append(j);
    }

    cb->E=cb0->E;
    ARRAY<IV3> E=cb1->E;
    E.Flattened()=index_v_map.Subset(E.Flattened());
    cb->E.Append_Elements(E);
    cb->Compute_Element_Edges();

    return *pr.x;
}
//#####################################################################
// Function Merge_Blocks
//#####################################################################
template<class T> template <class F> void COMPONENT_LAYOUT_FEM<T>::
Merge_Blocks(BLOCK_ID id,CON_ID con_id,BLOCK_ID id2,F alias)
{
    PHYSBAM_ASSERT(con_id>=CON_ID());
    PHYSBAM_ASSERT(id>=BLOCK_ID());
    BLOCK<T>& bl=blocks(id);
    CON_ID con_id2=bl.connections(con_id).con_id;
    BLOCK<T>& bl2=blocks(id2);

    auto& pr=Merge_Canonical_Blocks(bl.block,con_id,bl.xform,bl2.block,con_id2,
        bl2.xform);
    bl.block=pr.x;

    for(CON_ID c=con_id+1;c<bl.connections.m;c++)
    {
        const auto& cn=bl.connections(c);
        if(cn.is_regular) blocks(alias(cn.id)).connections(cn.con_id).con_id=c-1;
        else irregular_connections(cn.irreg_id).con_id=c-1;
        bl.connections(c-1)=bl.connections(c);
    }
    bl.connections.Pop();

    for(CON_ID c(0);c<bl2.connections.m;c++)
    {
        if(c==con_id2) continue;
        auto& cn=bl2.connections(c);
        if(cn.is_regular)
            blocks(alias(cn.id)).connections(cn.con_id).con_id=bl.connections.m;
        else irregular_connections(cn.irreg_id).con_id=bl.connections.m;
        bl.connections.Append(cn);
    }

    for(auto i:bl2.edge_on)
    {
        auto& ic=irregular_connections(i.x);
        auto& eo=ic.edge_on(i.y);
        eo.b=id;
        eo.e=pr.z(eo.e);
        eo.v0=pr.y(eo.v0);
        eo.v1=pr.y(eo.v1);
        bl.edge_on.Append(i);
    }

    bl2={0};
}
//#####################################################################
// Function Merge_Blocks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Merge_Blocks()
{
    UNION_FIND<BLOCK_ID> uf(blocks.m);
    HASHTABLE<BLOCK_ID,BLOCK_ID> ufs_root_blk;
    auto alias=[&uf,&ufs_root_blk](BLOCK_ID b)
    {
        return ufs_root_blk.Get_Default(uf.Find(b),b);
    };
    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        if(!blocks(b).block) continue;
        if(int mask=Separates_Dofs(b))
        {
            PHYSBAM_ASSERT(!(blocks(b).flags&1));
            CON_ID besti(-1);
            int best=INT_MAX;
            BLOCK_ID id2(-1);
            for(CON_ID i(0);i<blocks(b).connections.m;i++)
            {
                if(!(mask&(1<<Value(i)))) continue;
                if(!blocks(b).connections(i).is_regular) continue;
                BLOCK_ID d=alias(blocks(b).connections(i).id);
                PHYSBAM_ASSERT(blocks(d).block);
                if(blocks(d).flags&1) continue;
                int c=Approx_Dof_Count(d);
                if(c<best) // keep blocks small
                {
                    best=c;
                    besti=i;
                    id2=d;
                }
            }
            PHYSBAM_ASSERT(besti>=CON_ID());
            PHYSBAM_ASSERT(id2>=BLOCK_ID());
            Merge_Blocks(b,besti,id2,alias);
            uf.Union(b,id2);
            ufs_root_blk.Set(uf.Find(b),b);
            b--; // repeat the check on this block
        }
    }

    // Remap the block list to compress it
    BLOCK_ID next(0);
    ARRAY<BLOCK_ID,BLOCK_ID> mapping(blocks.m,use_init,BLOCK_ID(-11));
    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        BLOCK_ID p=uf.Find(b);
        if(mapping(p)<BLOCK_ID()) mapping(p)=next++;
        mapping(b)=mapping(p);
        if(blocks(b).block) blocks(mapping(b))=blocks(b);
    }

    blocks.Resize(next);
    for(auto& bl:blocks)
    {
        for(auto& c:bl.connections)
            if(c.is_regular)
            {
                assert(Value(mapping(c.id))>=0);
                c.id=mapping(c.id);
            }
    }
    for(auto& ic:irregular_connections)
    {
        ic.regular=mapping(ic.regular);
        for(auto& p:ic.edge_on)
            p.b=mapping(p.b);
    }
    for(auto& bc:bc_v) bc.b=mapping(bc.b);
    for(auto& bc:bc_t) bc.b=mapping(bc.b);
    assert(!reference_block_data.m);
    assert(!reference_connection_data.m);
    assert(!reference_irregular_data.m);
}
//#####################################################################
// Function Approx_Dof_Count
//#####################################################################
template<class T> int COMPONENT_LAYOUT_FEM<T>::
Approx_Dof_Count(BLOCK_ID b)
{
    const auto& bl=blocks(b);
    const auto* cb=bl.block;
    int num=cb->X.m;
    for(const auto& cs:cb->cross_sections)
        num-=cs.v.Size()/2;
    return num;
}

//#####################################################################
// Function Compute_Block_Hash
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Compute_Reference_Blocks()
{
    typedef ARRAY<std::tuple<CANONICAL_BLOCK<T>*,CON_ID,bool> > REG_CON;
    typedef ARRAY<std::tuple<CANONICAL_BLOCK<T>*,int> > IRREG_CON;
    typedef ARRAY<std::tuple<CANONICAL_BLOCK<T>*,CON_ID,int,int> > EDGE_ON;
    typedef std::tuple<CANONICAL_BLOCK<T>*,REG_CON,IRREG_CON> KEY;

    HASHTABLE<KEY,REFERENCE_BLOCK_ID> h;

    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        REG_CON reg_con;
        IRREG_CON irreg_con;
        EDGE_ON edge_on;

        auto& bl=blocks(b);
        for(CON_ID cc(0);cc<bl.connections.m;cc++)
        {
            const auto& c=bl.connections(cc);
            if(c.is_regular)
                reg_con.Append(std::make_tuple(blocks(c.id).block,c.con_id,c.master));
            else
            {
                reg_con.Append(std::make_tuple((CANONICAL_BLOCK<T>*)0,CON_ID(-1),true));

                const auto& ic=irregular_connections(c.irreg_id);
                Visit_Irregular_Cross_Section_Dofs(blocks,ic,
                    [&](const IRREGULAR_VISITOR& iv)
                    {
                        irreg_con.Append(std::make_tuple(blocks(iv.b).block,iv.ie));
                    });
            }
        }

        for(const auto p:bl.edge_on)
        {
            auto& ic=irregular_connections(p.x);
            edge_on.Append(std::make_tuple(blocks(ic.regular).block,ic.con_id,p.y,ic.edge_on(p.y).e));
        }

        auto pr=h.Insert(std::make_tuple(bl.block,reg_con,irreg_con),{});
        if(force_blk_ref || pr.y)
        {
            *pr.x=reference_block_data.Add_End();
            auto& rd=reference_block_data(*pr.x);
            rd.b=b;
        }
        bl.ref_id=*pr.x;
    }
}
//#####################################################################
// Function Compute_Reference_Regular_Connections
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Compute_Reference_Regular_Connections()
{
    for(BLOCK_ID b(0);b<blocks.m;b++)
    {
        const auto& bl=blocks(b);
        for(CON_ID cc(0);cc<bl.connections.m;cc++)
        {
            const auto& c=bl.connections(cc);
            if(c.is_regular && c.master)
            {
                const auto& other=blocks(c.id);
                REFERENCE_IRREGULAR_ID ic_ref[4]=
                {
                    bl.edge_on.m?irregular_connections(bl.edge_on(0).x).ref_id:REFERENCE_IRREGULAR_ID(-1),
                    bl.edge_on.m?irregular_connections(bl.edge_on(1).x).ref_id:REFERENCE_IRREGULAR_ID(-1),
                    other.edge_on.m?irregular_connections(other.edge_on(0).x).ref_id:REFERENCE_IRREGULAR_ID(-1),
                    other.edge_on.m?irregular_connections(other.edge_on(1).x).ref_id:REFERENCE_IRREGULAR_ID(-1)
                };
                auto pr=regular_connection_hash.Insert(std::make_tuple(
                    bl.ref_id,cc,ic_ref[0],ic_ref[1],
                    other.ref_id,c.con_id,ic_ref[2],ic_ref[3]),{});
                if(pr.y) *pr.x=reference_connection_data.Append({{b,c.id},{cc,c.con_id}});
            }
        }
    }
}
//#####################################################################
// Function Regular_Connection
//#####################################################################
template<class T> REFERENCE_CONNECTION_ID COMPONENT_LAYOUT_FEM<T>::
Regular_Connection(BLOCK_ID b0,CON_ID c0,BLOCK_ID b1) const
{
    const auto& bl=blocks(b0);
    const auto& other=blocks(b1);
    REFERENCE_IRREGULAR_ID ic_ref[4]=
    {
        bl.edge_on.m?irregular_connections(bl.edge_on(0).x).ref_id:REFERENCE_IRREGULAR_ID(-1),
        bl.edge_on.m?irregular_connections(bl.edge_on(1).x).ref_id:REFERENCE_IRREGULAR_ID(-1),
        other.edge_on.m?irregular_connections(other.edge_on(0).x).ref_id:REFERENCE_IRREGULAR_ID(-1),
        other.edge_on.m?irregular_connections(other.edge_on(1).x).ref_id:REFERENCE_IRREGULAR_ID(-1)
    };
    return regular_connection_hash.Get(std::make_tuple(
        bl.ref_id,c0,ic_ref[0],ic_ref[1],
        other.ref_id,bl.connections(c0).con_id,ic_ref[2],ic_ref[3]));
}
//#####################################################################
// Function Compute_Reference_Irregular_Connections
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Compute_Reference_Irregular_Connections()
{
    typedef std::tuple<CANONICAL_BLOCK<T>*,CON_ID,ARRAY<PAIR<CANONICAL_BLOCK<T>*,int> > > KEY;
    HASHTABLE<KEY,REFERENCE_IRREGULAR_ID> ref;

    for(IRREG_ID i(0);i<irregular_connections.m;i++)
    {
        auto& ic=irregular_connections(i);
        ARRAY<PAIR<CANONICAL_BLOCK<T>*,int> > ae;
        for(auto p:ic.edge_on) ae.Append({blocks(p.b).block,p.e});
        auto pr=ref.Insert(std::make_tuple(blocks(ic.regular).block,ic.con_id,ae),{});
        if(force_blk_ref || pr.y) *pr.x=reference_irregular_data.Append({i});
        ic.ref_id=*pr.x;
    }
}
//#####################################################################
// Function Compute_Connection_Hash
//#####################################################################
template<class T> int COMPONENT_LAYOUT_FEM<T>::
Compute_Connection_Hash(BLOCK_ID b0,CON_ID con_id0,BLOCK_ID b1,CON_ID con_id1)
{
    auto id0=blocks(b0).block,id1=blocks(b1).block;
    return Hash(std::make_tuple(id0,con_id0,id1,con_id1));
}
//#####################################################################
// Function Compute_Dof_Remapping
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Compute_Dof_Remapping(REFERENCE_BLOCK_DATA& rd)
{
    const auto& bl=blocks(rd.b);
    const auto* cb=bl.block;
    rd.dof_map_v.Resize(cb->X.m,use_init,1);
    rd.dof_map_e.Resize(cb->S.m,use_init,1);
    
    for(auto i:bl.edge_on)
    {
        auto& ic=irregular_connections(i.x);
        auto& id=ic.edge_on(i.y);
        if(i.y<ic.edge_on.m/2+1)
            rd.dof_map_v(id.v0)=0;
        if(i.y<ic.edge_on.m/2)
        {
            rd.dof_map_v(id.v1)=0;
            rd.dof_map_e(id.e)=0;
        }
    }

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
        {
            assert(blocks(c.id).connections(c.con_id).id==rd.b);
            assert(blocks(c.id).connections(c.con_id).con_id==cc);
            assert(blocks(c.id).connections(c.con_id).master!=c.master);
            const auto* cb2=blocks(c.id).block;
            Visit_Regular_Cross_Section_Dofs(cb->cross_sections(cc),
                cb2->cross_sections(c.con_id),c.master,
                [&](int a,int b,bool o)
                {
                    if(!o) rd.dof_map_v(a)=0;
                },
                [&](int a,int b,bool o)
                {
                    if(!o) rd.dof_map_e(a)=0;
                });
        }
        else
        {
            const auto& ic=irregular_connections(c.irreg_id);
            Visit_Irregular_Cross_Section_Dofs(blocks,ic,
                [&](const IRREGULAR_VISITOR& iv)
                {
                    if(!iv.be) rd.dof_map_e(iv.re)=0;
                    if(!iv.b0) rd.dof_map_v(iv.r0)=0;
                    if(!iv.b1) rd.dof_map_v(iv.r1)=0;
                });
        }
    }

    int iv=0,ie=0,ip=0;
    rd.dof_map_p=rd.dof_map_v;
    rd.dof_map_v.Subset(cb->bc_v).Fill(0);
    rd.dof_map_e.Subset(cb->bc_e).Fill(0);
    for(int i=0;i<rd.dof_map_v.m;i++)
        rd.dof_map_v(i)=rd.dof_map_v(i)?iv++:-1;
    for(int i=0;i<rd.dof_map_e.m;i++)
        rd.dof_map_e(i)=rd.dof_map_e(i)?ie++:-1;
    for(int i=0;i<rd.dof_map_p.m;i++)
        rd.dof_map_p(i)=rd.dof_map_p(i)?ip++:-1;
    rd.num_dofs_s={cb->X.m,cb->S.m,cb->X.m};
    rd.num_dofs_d={iv,ie,ip};
}
//#####################################################################
// Function Compute_Bounding_Box
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<T>::
Compute_Bounding_Box() const -> RANGE<TV>
{
    RANGE<TV> box=RANGE<TV>::Empty_Box();
    for(const auto&bl:blocks)
        for(auto i:bl.block->bc_v)
            box.Enlarge_To_Include_Point(bl.xform*bl.block->X(i));
    return box;
}
//#####################################################################
// Function Compute_Dof_Pairs
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Compute_Dof_Pairs(REFERENCE_BLOCK_DATA& rd)
{
    // Copy from self
    for(int i=0;i<rd.dof_map_v.m;i++)
        if(rd.dof_map_v(i)>=0)
            rd.pairs.v.Append({rd.dof_map_v(i),i});
    for(int i=0;i<rd.dof_map_e.m;i++)
        if(rd.dof_map_e(i)>=0)
            rd.pairs.e.Append({rd.dof_map_e(i),i});
    for(int i=0;i<rd.dof_map_p.m;i++)
        if(rd.dof_map_p(i)>=0)
            rd.pairs.p.Append({rd.dof_map_p(i),i});
    Fill_Num_Dofs(rd.pairs,rd.b,rd.b);
}
//#####################################################################
// Function Compute_Dof_Pairs
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Compute_Dof_Pairs(REFERENCE_CONNECTION_DATA& rc)
{
    const auto* cb0=blocks(rc.b[0]).block;
    const auto* cb1=blocks(rc.b[1]).block;
    REFERENCE_BLOCK_DATA* rd[2];
    rd[0]=&reference_block_data(blocks(rc.b[0]).ref_id);
    rd[1]=&reference_block_data(blocks(rc.b[1]).ref_id);
    Fill_Num_Dofs(rc.reg_pairs[0],rc.b[0],rc.b[1]);
    Fill_Num_Dofs(rc.reg_pairs[1],rc.b[1],rc.b[0]);
    Visit_Regular_Cross_Section_Dofs(cb0->cross_sections(rc.con_id[0]),
        cb1->cross_sections(rc.con_id[1]),true, // 0 is master
        [&](int a,int b,bool o)
        {
            int f[2]={a,b};
            int s=o,d=!o,v=rd[d]->dof_map_v(f[d]),p=rd[d]->dof_map_p(f[d]);
            if(v>=0) rc.reg_pairs[d].v.Append({v,f[s]});
            if(p>=0) rc.reg_pairs[d].p.Append({p,f[s]});
        },
        [&](int a,int b,bool o)
        {
            int f[2]={a,b};
            int s=o,d=!o,e=rd[d]->dof_map_e(f[d]);
            if(e>=0) rc.reg_pairs[d].e.Append({e,f[s]});
        });
}
//#####################################################################
// Function Compute_Dof_Pairs
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Compute_Dof_Pairs(REFERENCE_IRREGULAR_DATA& ri)
{
    auto& ic=irregular_connections(ri.ic_id);

    // Copy from irregular neighbors
    REFERENCE_BLOCK_DATA* rd[2];
    rd[0]=&reference_block_data(blocks(ic.regular).ref_id);
    RID_ID id(-1);
    Visit_Irregular_Cross_Section_Dofs(blocks,ic,
        [&](const IRREGULAR_VISITOR& iv)
        {
            if(iv.n0) id=ri.pairs.Add_End();
            ri.mapping.Append({id,iv.n0});
            rd[1]=&reference_block_data(blocks(iv.b).ref_id);
            auto& z=ri.pairs(id);
            z.b=iv.b;

            {
                int s=iv.be,d=!s;
                DOF_PAIRS& q=z.irreg_pairs[d];
                int f[2]={iv.re,iv.ie};
                int r=rd[d]->dof_map_e(f[d]);
                if(r>=0) q.e.Append({r,f[s]});
            }

            if(iv.n0)
            {
                int s=iv.b0,d=!s;
                DOF_PAIRS& q=z.irreg_pairs[d];
                int f[2]={iv.r0,iv.i0};
                int r=rd[d]->dof_map_v(f[d]),p=rd[d]->dof_map_p(f[d]);
                if(r>=0) q.v.Append({r,f[s]});
                if(p>=0) q.p.Append({p,f[s]});
            }

            {
                int s=iv.b1,d=!s;
                DOF_PAIRS& q=z.irreg_pairs[d];
                int f[2]={iv.r1,iv.i1};
                int r=rd[d]->dof_map_v(f[d]),p=rd[d]->dof_map_p(f[d]);
                if(r>=0) q.v.Append({r,f[s]});
                if(p>=0) q.p.Append({p,f[s]});
            }
        });

    for(auto& r:ri.pairs)
    {
        Fill_Num_Dofs(r.irreg_pairs[0],ic.regular,r.b);
        Fill_Num_Dofs(r.irreg_pairs[1],r.b,ic.regular);
    }
}
//#####################################################################
// Function Compute_Dof_Pairs
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Compute_Dof_Pairs()
{
    Compute_Reference_Blocks();
    Compute_Reference_Irregular_Connections();
    Compute_Reference_Regular_Connections();

    for(auto& rd:reference_block_data)
    {
        Compute_Dof_Remapping(rd);
        Compute_Dof_Pairs(rd);
    }

    for(auto& rc:reference_connection_data)
        Compute_Dof_Pairs(rc);

    for(auto& ri:reference_irregular_data)
        Compute_Dof_Pairs(ri);
}
//#####################################################################
// Function Fill_Num_Dofs
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Fill_Num_Dofs(DOF_PAIRS& dp,BLOCK_ID d,BLOCK_ID s)
{
    dp.num_dofs_d=reference_block_data(blocks(d).ref_id).num_dofs_d;
    dp.num_dofs_s=reference_block_data(blocks(s).ref_id).num_dofs_s;
}
//#####################################################################
// Function Regular_Connection_Pair
//#####################################################################
template<class T> auto COMPONENT_LAYOUT_FEM<T>::
Regular_Connection_Pair(BLOCK_ID b,CON_ID con_id,bool is_dest) -> const DOF_PAIRS&
{
    auto& bl=blocks(b);
    BLOCK_ID b1=bl.connections(con_id).id;
    CON_ID con_id1=bl.connections(con_id).con_id;

    if(bl.connections(con_id).master)
    {
        auto ref_con_id=Regular_Connection(b,con_id,b1);
        return reference_connection_data(ref_con_id).reg_pairs[!is_dest];
    }
    auto ref_con_id=Regular_Connection(b1,con_id1,b);
    return reference_connection_data(ref_con_id).reg_pairs[is_dest];
}
//#####################################################################
// Function Fill_Reference_Ticks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Fill_Reference_Ticks()
{
    HASHTABLE<const CANONICAL_BLOCK<T>*,ARRAY<bool> > half_edge_dir; // ccw->false
    for(const auto& bl:blocks)
    {
        const auto* cb=bl.block;
        auto pr=half_edge_dir.Insert(cb,{});
        if(!pr.y) continue;
        auto& ee=*pr.x;
        ee.Resize(cb->S.m);
        for(int i=0;i<cb->E.m;i++)
        {
            for(int j=0;j<3;j++)
            {
                const auto& edge=cb->element_edges(i)(j);
                ee(edge.x)=edge.y;
            }
        }
    }

    for(auto& rd:reference_block_data)
        rd.ticks_e=blocks(rd.b).block->ticks;

    for(auto& rd:reference_block_data)
    {
        const auto& bl=blocks(rd.b);
        const auto* cb=bl.block;
        for(CON_ID cc(0);cc<bl.connections.m;cc++)
        {
            const auto& c=bl.connections(cc);
            if(c.is_regular)
            {
                const auto* cb2=blocks(c.id).block;
                Visit_Regular_Cross_Section_Dofs(cb->cross_sections(cc),
                    cb2->cross_sections(c.con_id),c.master,
                    [&](int a,int b,bool o){},
                    [&](int a,int b,bool o)
                    {
                        int ta=rd.ticks_e(a);
                        int& tb=reference_block_data(blocks(c.id).ref_id).ticks_e(b);
                        bool ccw_a=half_edge_dir.Get(cb)(a),ccw_b=half_edge_dir.Get(cb2)(b);
                        if((ccw_a!=ccw_b && ta!=tb) || (ccw_a==ccw_b && ta==tb))
                            tb=1-tb;
                    });
            }
            else
            {
                const auto& ic=irregular_connections(c.irreg_id);
                Visit_Irregular_Cross_Section_Dofs(blocks,ic,
                    [&](const IRREGULAR_VISITOR& iv)
                    {
                        int& t_reg=rd.ticks_e(iv.re);
                        int t_irreg=reference_block_data(blocks(iv.b).ref_id).ticks_e(iv.ie);
                        bool ccw_reg=half_edge_dir.Get(cb)(iv.re),ccw_irreg=half_edge_dir.Get(blocks(iv.b).block)(iv.ie);
                        if((ccw_reg!=ccw_irreg && t_reg!=t_irreg) || (ccw_reg==ccw_irreg && t_reg==t_irreg))
                            t_reg=1-t_reg;
                    });
            }
        }
    }

    Fill_Element_Tick_Masks();
}
//#####################################################################
// Function Fill_Element_Tick_Masks
//#####################################################################
template<class T> void COMPONENT_LAYOUT_FEM<T>::
Fill_Element_Tick_Masks()
{
    for(auto& rb:reference_block_data)
    {
        rb.ticks_t.Resize(blocks(rb.b).block->E.m);
        const auto* cb=blocks(rb.b).block;
        for(int i=0;i<cb->E.m;i++)
        {
            int m=0;
            for(int j=0;j<3;j++)
            {
                auto p=cb->element_edges(i)(j);
                m|=(rb.ticks_e(p.x)^p.y)<<j;
            }
            rb.ticks_t(i)=m;
        }
    }
}
template class COMPONENT_LAYOUT_FEM<double>;
}
