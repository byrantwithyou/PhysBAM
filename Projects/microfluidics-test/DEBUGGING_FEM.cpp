//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Images/EPS_FILE.h>
#include "DEBUGGING_FEM.h"
#include "MATRIX_CONSTRUCTION_FEM.h"
#include "VISITORS_FEM.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> DEBUGGING_FEM<T>::
DEBUGGING_FEM(COMPONENT_LAYOUT_FEM<T>& cl)
    :cl(cl)
{
}
//#####################################################################
// Function Visualize_Block_State
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Block_State(BLOCK_ID b) const
{
    auto& bl=cl.blocks(b);
    auto* cb=bl.block;
    auto Z=[=](int i){return bl.xform*cb->X(i);};
    for(auto t:cb->E)
    {
        VECTOR<TV2,3> P;
        for(int i=0;i<3;i++) P(i)=bl.xform*cb->X(t(i));
        for(auto p:P) Add_Debug_Object(VECTOR<TV2,2>(p,P.Average()),VECTOR<T,3>(.5,.5,.5));
    }
    HASHTABLE<int> he,hv;
    he.Set_All(cb->bc_e);
    hv.Set_All(cb->bc_v);
    for(int i=0;i<cb->S.m;i++)
    {
        TV2 X=Z(cb->S(i).x);
        TV2 Y=Z(cb->S(i).y);
        Add_Debug_Object(VECTOR<TV2,2>(X,Y),he.Contains(i)?VECTOR<T,3>(0,1,1):VECTOR<T,3>(0,0,1));
    }
    for(int i=0;i<cb->X.m;i++)
    {
        Add_Debug_Particle(Z(i),hv.Contains(i)?VECTOR<T,3>(1,1,0):VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV2>("display_size",.1);
    }

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
        {
            const auto* cb2=cl.blocks(c.id).block;
            Visit_Regular_Cross_Section_Dofs(cb->cross_sections(cc),
                cb2->cross_sections(c.con_id),c.master,
                [&](int a,int b,bool o)
                {
                    Add_Debug_Particle(Z(a),o?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.2);
                },
                [&](int a,int b,bool o)
                {
                    Add_Debug_Particle((Z(cb->S(a).x)+Z(cb->S(a).y))/2,o?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.2);
                });
        }
        else
        {
            const auto& ic=cl.irregular_connections(c.irreg_id);
            Visit_Irregular_Cross_Section_Dofs(cl.blocks,ic,
                [&](const IRREGULAR_VISITOR& iv)
                {
                    TV2 A=Z(iv.r0),B=Z(iv.r1);
                    Add_Debug_Particle(A,VECTOR<T,3>(iv.b0,0,1));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.2);
                    Add_Debug_Particle(B,VECTOR<T,3>(iv.b1,0,1));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.25);
                    Add_Debug_Particle((A+B)/2,VECTOR<T,3>(iv.be,0,1));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.2);
                });
        }
    }

    for(auto i:bl.edge_on)
    {
        auto& ic=cl.irregular_connections(i.x);
        auto& id=ic.edge_on(i.y);
        bool o0=i.y>=(ic.edge_on.m/2+1);
        bool o1=(i.y>=ic.edge_on.m/2);

        TV2 A=Z(id.v0),B=Z(id.v1);
        Add_Debug_Particle(A,VECTOR<T,3>(1,1,1)/(1+o0));
        Debug_Particle_Set_Attribute<TV2>("display_size",.05);
        Add_Debug_Particle(B,VECTOR<T,3>(1,1,1)/(1+o1));
        Debug_Particle_Set_Attribute<TV2>("display_size",.06);
        Add_Debug_Particle((A+B)/2,VECTOR<T,3>(1,1,1)/(1+o1));
        Debug_Particle_Set_Attribute<TV2>("display_size",.05);
    }
}
//#####################################################################
// Function Visualize_Block_State
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Block_Dofs(BLOCK_ID b) const
{
    auto& bl=cl.blocks(b);
    auto* cb=bl.block;
    const auto& rd=cl.reference_block_data(cl.blocks(b).ref_id);
    auto Z=[=](int i){return bl.xform*cb->X(i);};
    for(auto t:cb->E)
    {
        VECTOR<TV2,3> P;
        for(int i=0;i<3;i++) P(i)=bl.xform*cb->X(t(i));
        for(auto p:P) Add_Debug_Object(VECTOR<TV2,2>(p,P.Average()),VECTOR<T,3>(.5,.5,.5));
    }
    HASHTABLE<int> he,hv;
    he.Set_All(cb->bc_e);
    hv.Set_All(cb->bc_v);
    for(int i=0;i<cb->S.m;i++)
    {
        TV2 X=Z(cb->S(i).x);
        TV2 Y=Z(cb->S(i).y);
        Add_Debug_Object(VECTOR<TV2,2>(X,Y),he.Contains(i)?VECTOR<T,3>(0,1,1):VECTOR<T,3>(0,0,1));
    }
    Flush_Frame<TV2>(LOG::sprintf("block %P\n",b).c_str());
    for(int i=0;i<cb->X.m;i++)
    {
        Add_Debug_Text(Z(i),LOG::sprintf("%d",i),VECTOR<T,3>(1,1,0));
    }
    for(int i=0;i<cb->S.m;i++)
    {
        Add_Debug_Text((Z(cb->S(i).x)+Z(cb->S(i).y))/2,LOG::sprintf("%d",i),VECTOR<T,3>(0,1,0));
    }
    Flush_Frame<TV2>(LOG::sprintf("block %P full\n",b).c_str());
    for(int i=0;i<cb->X.m;i++)
    {
        int k=rd.dof_map_v(i);
        if(k>=0) Add_Debug_Text(Z(i),LOG::sprintf("%d",k),VECTOR<T,3>(1,1,0));
    }
    Flush_Frame<TV2>(LOG::sprintf("block %P v\n",b).c_str());
    for(int i=0;i<cb->S.m;i++)
    {
        int k=rd.dof_map_e(i);
        if(k>=0) Add_Debug_Text((Z(cb->S(i).x)+Z(cb->S(i).y))/2,LOG::sprintf("%d",k),VECTOR<T,3>(0,1,0));
    }
    Flush_Frame<TV2>(LOG::sprintf("block %P e\n",b).c_str());
    for(int i=0;i<cb->X.m;i++)
    {
        int k=rd.dof_map_p(i);
        if(k>=0) Add_Debug_Text(Z(i),LOG::sprintf("%d",k),VECTOR<T,3>(1,1,0));
    }
    Flush_Frame<TV2>(LOG::sprintf("block %P p\n",b).c_str());

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
        {
            const auto* cb2=cl.blocks(c.id).block;
            Visit_Regular_Cross_Section_Dofs(cb->cross_sections(cc),
                cb2->cross_sections(c.con_id),c.master,
                [&](int a,int b,bool o)
                {
                    Add_Debug_Particle(Z(a),o?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.2);
                },
                [&](int a,int b,bool o)
                {
                    Add_Debug_Particle((Z(cb->S(a).x)+Z(cb->S(a).y))/2,o?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.2);
                });
        }
        else
        {
            const auto& ic=cl.irregular_connections(c.irreg_id);
            Visit_Irregular_Cross_Section_Dofs(cl.blocks,ic,
                [&](const IRREGULAR_VISITOR& iv)
                {
                    TV2 A=Z(iv.r0),B=Z(iv.r1);
                    Add_Debug_Particle(A,VECTOR<T,3>(iv.b0,0,1));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.2);
                    Add_Debug_Particle(B,VECTOR<T,3>(iv.b1,0,1));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.25);
                    Add_Debug_Particle((A+B)/2,VECTOR<T,3>(iv.be,0,1));
                    Debug_Particle_Set_Attribute<TV2>("display_size",.2);
                });
        }
    }

    for(auto i:bl.edge_on)
    {
        auto& ic=cl.irregular_connections(i.x);
        auto& id=ic.edge_on(i.y);
        bool o0=i.y>=(ic.edge_on.m/2+1);
        bool o1=(i.y>=ic.edge_on.m/2);

        TV2 A=Z(id.v0),B=Z(id.v1);
        Add_Debug_Particle(A,VECTOR<T,3>(1,1,1)/(1+o0));
        Debug_Particle_Set_Attribute<TV2>("display_size",.05);
        Add_Debug_Particle(B,VECTOR<T,3>(1,1,1)/(1+o1));
        Debug_Particle_Set_Attribute<TV2>("display_size",.06);
        Add_Debug_Particle((A+B)/2,VECTOR<T,3>(1,1,1)/(1+o1));
        Debug_Particle_Set_Attribute<TV2>("display_size",.05);
    }
}
//#####################################################################
// Function Visualize_Solution
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Solution(const BLOCK_VECTOR<TV3>& U,BLOCK_ID b,bool remap_dofs) const
{
    if(!remap_dofs) return;

    const auto& rb=cl.reference_block_data(cl.blocks(b).ref_id);
    const auto& M=cl.blocks(b).xform;
    DOF_LAYOUT<TV3> dl(cl,rb,true);

    Visualize_Tetrahedron(b);
    Visit_Compressed_Dofs(dl,rb,
        [&M,&U](int v,const TV3& X)
        {
            Add_Debug_Particle(xform(M,X),VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV3>("V",U.Get_v(v));
        },
        [&M,&U](int e,const TV3& X)
        {
            Add_Debug_Particle(xform(M,X),VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV3>("V",U.Get_e(e));
        },
        [](int p,const TV3& X){});
    Flush_Frame<TV3>(LOG::sprintf("sol %P ve",b).c_str());

    Visualize_Tetrahedron(b);
    Visit_Compressed_Dofs(dl,rb,
        [](int v,const TV3& X){},
        [](int e,const TV3& X){},
        [&M,&U](int p,const TV3& X)
        {
            Add_Debug_Particle(xform(M,X),VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV3>("display_size",U.Get_p(p));
        });
    Flush_Frame<TV3>(LOG::sprintf("sol %P p",b).c_str());
}
//#####################################################################
// Function Visualize_Solution
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Solution(const BLOCK_VECTOR<TV2>& U,BLOCK_ID b,bool remap_dofs) const
{
    const auto& bl=cl.blocks(b);
    const auto* cb=bl.block;
    auto Z=[=](int i){return bl.xform*cb->X(i);};
    const auto& rd=cl.reference_block_data(cl.blocks(b).ref_id);
    for(int i=0;i<cb->X.m;i++)
    {
        int k=remap_dofs?rd.dof_map_v(i):i;
        if(k>=0)
        {
            Add_Debug_Particle(Z(i),VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV2>("V",U.Get_v(k));
        }
    }

    for(int i=0;i<cb->S.m;i++)
    {
        int k=remap_dofs?rd.dof_map_e(i):i;
        if(k>=0)
        {
            Add_Debug_Particle((Z(cb->S(i).x)+Z(cb->S(i).y))/2,VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV2>("V",U.Get_e(k));
        }
    }

    for(int i=0;i<cb->X.m;i++)
    {
        int k=remap_dofs?rd.dof_map_p(i):i;
        if(k>=0)
        {
            Add_Debug_Particle(Z(i),VECTOR<T,3>(0,1,0));
            Debug_Particle_Set_Attribute<TV2>("display_size",U.Get_p(k));
        }
    }
}
//#####################################################################
// Function Visualize_Ticks
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Ticks(BLOCK_ID b,bool reference_ticks) const
{
    auto& bl=cl.blocks(b);
    auto* cb=bl.block;
    auto Z=[=](int i){return bl.xform*cb->X(i);};

    HASHTABLE<IV2,int> edge_lookup;
    for(int i=0;i<cb->S.m;i++)
        edge_lookup.Set(cb->S(i).Sorted(),i);

    const ARRAY<int>& ticks=reference_ticks?cl.reference_block_data(bl.ref_id).ticks_e:cb->ticks;
    for(int j=0;j<cb->E.m;j++)
    {
        const auto& t=cb->E(j);
        int check[3]={0,0,0};
        for(int i=0;i<3;i++)
        {
            int e0=cb->element_edges(j)((i+1)%3).x,e1=cb->element_edges(j)((i+2)%3).x;
            TV2 X=Z(t(i)),X0=Z(t((i+2)%3)),X1=Z(t((i+1)%3));
            TV2 x0=0.85*X+0.13*X0+0.02*X1,x1=0.85*X+0.02*X0+0.13*X1;
            int cnt=0;
            if(cb->S(e0)(ticks(e0))==t(i))
            {
                cnt++;
                Add_Debug_Particle(x0,VECTOR<T,3>(1,0,0));
            }
            if(cb->S(e1)(ticks(e1))==t(i))
            {
                cnt++;
                Add_Debug_Particle(x1,VECTOR<T,3>(1,0,0));
            }
            check[cnt]++;
        }
        PHYSBAM_ASSERT(check[0]==1 && check[1]==1 && check[2]==1);
    }
}
//#####################################################################
// Function Visualize_Flat_Dofs
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Flat_Dofs() const
{
    int next_u=0,next_p=0;
    ARRAY<int,BLOCK_ID> first[3];
    for(int i=0;i<3;i++) first[i].Resize(cl.blocks.m);
    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& c=cl.reference_block_data(cl.blocks(b).ref_id);
        first[0](b)=next_u;
        next_u+=c.num_dofs_d.v*2;
        first[1](b)=next_u;
        next_u+=c.num_dofs_d.e*2;
        first[2](b)=next_p;
        next_p+=c.num_dofs_d.p;
    }
    first[2]+=next_u;

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        Visualize_Block_State(b);
        auto& bl=cl.blocks(b);
        const auto* cb=bl.block;
        auto Z=[=](int i){return bl.xform*cb->X(i);};
        const auto& rd=cl.reference_block_data(bl.ref_id);
        for(int i=0;i<cb->X.m;i++)
            if(rd.dof_map_v(i)>=0)
            {
                int dof=first[0](b)+rd.dof_map_v(i)*2;
                Add_Debug_Text(Z(i),LOG::sprintf("%d",dof),VECTOR<T,3>(1,1,0));
            }
        for(int i=0;i<cb->S.m;i++)
            if(rd.dof_map_e(i)>=0)
            {
                int dof=first[1](b)+rd.dof_map_e(i)*2;
                Add_Debug_Text((Z(cb->S(i).x)+Z(cb->S(i).y))/2,LOG::sprintf("%d",dof),VECTOR<T,3>(1,1,0));
            }
    }
    Flush_Frame<TV2>("flat velocity dofs");
    
    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        Visualize_Block_State(b);
        auto& bl=cl.blocks(b);
        const auto* cb=bl.block;
        auto Z=[=](int i){return bl.xform*cb->X(i);};
        const auto& rd=cl.reference_block_data(bl.ref_id);
        for(int i=0;i<cb->X.m;i++)
            if(rd.dof_map_p(i)>=0)
            {
                int dof=first[2](b)+rd.dof_map_p(i);
                Add_Debug_Text(Z(i),LOG::sprintf("%d",dof),VECTOR<T,3>(1,1,0));
            }
    }
    Flush_Frame<TV2>("flat pressure dofs");
    

}
//#####################################################################
// Function Visualize_Tetrahedron
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Tetrahedron(BLOCK_ID b) const
{
    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto* cb=cl.blocks(b).block;
        auto Z=[=](int i){return cl.blocks(b).xform*cb->X(i);};
        for(auto t:cb->E)
        {
            VECTOR<TV3,3> P;
            for(int i=0;i<3;i++) P(i)=Z(t(i)).Append(-0.01);
            for(int i=0;i<3;i++) Add_Debug_Object(VECTOR<TV3,2>(P(i),P((i+1)%3)),VECTOR<T,3>(.2,.2,.2));
        }
    }

    const auto& rb=cl.reference_block_data(cl.blocks(b).ref_id);
    const auto& M=cl.blocks(b).xform;
    DOF_LAYOUT<TV3> dl(cl,rb,false);
    Visit_Elements(dl,[&M](const VISIT_ELEMENT_DATA<TV3>& vt)
    {
        VECTOR<TV3,4> P;
        for(int i=0;i<4;i++) P(i)=xform(M,vt.X(i));
        Add_Debug_Object(VECTOR<TV3,2>(P(0),P(1)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV3,2>(P(0),P(2)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV3,2>(P(0),P(3)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV3,2>(P(1),P(2)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV3,2>(P(2),P(3)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV3,2>(P(3),P(1)),VECTOR<T,3>(1,1,1));
    });
}
//#####################################################################
// Function Visualize_Tetrahedron_Dofs
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Tetrahedron_Dofs(const MATRIX_CONSTRUCTION_FEM<TV3>& mc) const
{
    ARRAY<int,BLOCK_ID> first[3];
    mc.Compute_Global_Dof_Mapping(first);

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& rb=cl.reference_block_data(cl.blocks(b).ref_id);
        const auto& M=cl.blocks(b).xform;
        DOF_LAYOUT<TV3> dl(cl,rb,true);
        int count_v=0,count_e=0,count_p=0;

        Visualize_Tetrahedron(b);
        Visit_Compressed_Dofs(dl,rb,
            [&M,b,&count_v,&first](int v,const TV3& X)
            {
                Add_Debug_Particle(xform(M,X),VECTOR<T,3>(1,1,1));
                int dof=first[0](b)+count_v*TV3::m;
                Debug_Particle_Set_Attribute<TV3>("v",VECTOR<int,2>(v,dof));
                ++count_v;
            },
            [&M,b,&count_e,&first](int e,const TV3& X)
            {
                Add_Debug_Particle(xform(M,X),VECTOR<T,3>(1,1,1));
                int dof=first[1](b)+count_e*TV3::m;
                Debug_Particle_Set_Attribute<TV3>("e",VECTOR<int,2>(e,dof));
                ++count_e;
            },
            [](int p,const TV3& X){});
        Flush_Frame<TV3>(LOG::sprintf("block %P ve",b).c_str());

        Visualize_Tetrahedron(b);
        Visit_Compressed_Dofs(dl,rb,
            [](int v,const TV3& X){},
            [](int e,const TV3& X){},
            [&M,b,&count_p,&first](int p,const TV3& X)
            {
                Add_Debug_Particle(xform(M,X),VECTOR<T,3>(1,1,1));
                int dof=first[2](b)+count_p;
                Debug_Particle_Set_Attribute<TV3>("p",VECTOR<int,2>(p,dof));
                ++count_p;
            });
        Flush_Frame<TV3>(LOG::sprintf("block %P p",b).c_str());
    }
}
//#####################################################################
// Function Highlight_Dof
//#####################################################################
template<class T> template<int d> void DEBUGGING_FEM<T>::
Highlight_Dof(BLOCK_ID b,int vep,int r,int dim) const
{
    typedef VECTOR<T,d> TV;
    const auto& M=cl.blocks(b).xform;
    const auto& rb=cl.reference_block_data(cl.blocks(b).ref_id);
    DOF_LAYOUT<TV> dl(cl,rb,true);
    int count=0;
    VECTOR<T,3> color;
    if(dim>=0) color(dim)=1;
    Visit_Compressed_Dofs(dl,rb,
        [&count,&M,&color,vep,r](int v,const TV& X)
        {
            if(vep!=0) return;
            if(count++==r)
            {
                Add_Debug_Particle(xform(M,X),color);
                Debug_Particle_Set_Attribute<TV>("display_size",.5);
            }
        },
        [&count,&M,&color,vep,r](int e,const TV& X)
        {
            if(vep!=1) return;
            if(count++==r)
            {
                Add_Debug_Particle(xform(M,X),color);
                Debug_Particle_Set_Attribute<TV>("display_size",.5);
            }
        },
        [&count,&M,&color,vep,r](int p,const TV& X)
        {
            if(vep!=2) return;
            if(count++==r)
            {
                Add_Debug_Particle(xform(M,X),VECTOR<T,3>(1,1,0));
                Debug_Particle_Set_Attribute<TV>("display_size",.5);
            }
        });
}
//#####################################################################
// Function Visualize_Meshing
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Meshing(const std::string& name,const IV2& dim,const RANGE<TV2>& domain) const
{
    EPS_FILE<T> eps(name,RANGE<TV2>(TV2(),TV2(dim(0),dim(1))));
    T scale=(dim/domain.Edge_Lengths()).Min();
    eps.cur_format.line_width=1/scale;
    eps.cur_format.point_radius=1/scale;
    eps.cur_format.line_style=1;
    eps.cur_format.fill_style=1;

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& bl=cl.blocks(b);
        const auto* cb=bl.block;
        auto Z=[=](int i){return bl.xform*cb->X(i);};
        for(int i=0;i<cb->S.m;i++)
        {
            TV2 X=Z(cb->S(i).x);
            TV2 Y=Z(cb->S(i).y);
            if(!domain.Lazy_Inside(X) && !domain.Lazy_Inside(Y))
                continue;
            eps.Draw_Object(domain.Clamp(X)*scale,domain.Clamp(Y)*scale);
        }
    }
}
//#####################################################################
// Function Visualize_Domain
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Domain(const std::string& name,bool fill,const IV2& dim,const RANGE<TV2>& anno,const RANGE<TV2>& bound) const
{
    EPS_FILE<T> eps(name,RANGE<TV2>(TV2(),TV2(dim(0),dim(1))));
    RANGE<TV2> box=cl.Compute_Bounding_Box();
    T scale=(dim/box.Edge_Lengths()).Min();
    eps.cur_format.line_width=1/scale;
    eps.cur_format.point_radius=1/scale;

    if(bound.Empty())
    {
        for(const auto& bc:cl.bc_v)
        {
            eps.cur_format.line_style=0;
            eps.cur_format.fill_color=TV3(0,0,1);
            eps.cur_format.fill_style=1;
            const auto& bl=cl.blocks(bc.b);
            if(bc.flow_rate==0) continue;
            TV2 X;
            for(int i:bc.bc_v)
                X+=bl.xform*bl.block->X(i);
            X/=bc.bc_v.Size();
            eps.Draw_Object(X+bc.normal*12/scale,6/scale);
        }
        for(const auto& bc:cl.bc_t)
        {
            eps.cur_format.line_style=0;
            eps.cur_format.fill_color=TV3(1,0,0);
            eps.cur_format.fill_style=1;
            const auto& bl=cl.blocks(bc.b);
            TV2 X;
            for(int i:bc.bc_v)
                X+=bl.xform*bl.block->X(i);
            X/=bc.bc_v.Size();
            eps.Draw_Object(X+bc.normal*12/scale,6/scale);
        }
    }

    if(fill)
    {
        eps.cur_format.line_style=0;
        eps.cur_format.fill_style=1;
        eps.cur_format.fill_color=TV3(0,0,0);
        for(const auto& bl:cl.blocks)
            for(const auto& e:bl.block->E)
            {
                bool inside=true;
                for(int j:e)
                    if(!bound.Empty() && !bound.Lazy_Inside(bl.xform*bl.block->X(j))) inside=false;
                if(inside)
                    eps.Draw_Object(bl.xform*bl.block->X(e(0)),bl.xform*bl.block->X(e(1)),bl.xform*bl.block->X(e(2)));
            }
    }
    else
    {
        eps.cur_format.line_style=1;
        eps.cur_format.fill_style=0;
        auto draw_boundary_edge=[&eps,&bound](const BLOCK<T>& bl,int i)
        {
            const auto& e=bl.block->S(i);
            bool inside=true;
            for(int j:e)
                if(!bound.Empty() && !bound.Lazy_Inside(bl.xform*bl.block->X(j))) inside=false;
            if(inside)
            eps.Draw_Object(bl.xform*bl.block->X(e(0)),bl.xform*bl.block->X(e(1)));
        };

        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            const auto& bl=cl.blocks(b);
            for(int i:bl.block->bc_e)
                draw_boundary_edge(bl,i);
        }
        for(const auto& bc:cl.bc_v)
        {
            const auto& bl=cl.blocks(bc.b);
            for(int i:bc.bc_e)
                draw_boundary_edge(bl,i);
        }
        for(const auto& bc:cl.bc_t)
        {
            const auto& bl=cl.blocks(bc.b);
            for(int i:bc.bc_e)
                draw_boundary_edge(bl,i);
        }
    }

    if(!anno.Empty())
    {
        eps.cur_format.line_width=2/scale;
        eps.cur_format.line_style=1;
        eps.cur_format.line_color=TV3(1,0,0);
        eps.cur_format.fill_style=0;
        eps.Draw_Object(anno.min_corner,TV2(anno.max_corner.x,anno.min_corner.y));
        eps.Draw_Object(anno.min_corner,TV2(anno.min_corner.x,anno.max_corner.y));
        eps.Draw_Object(anno.max_corner,TV2(anno.max_corner.x,anno.min_corner.y));
        eps.Draw_Object(anno.max_corner,TV2(anno.min_corner.x,anno.max_corner.y));
    }
}
template class DEBUGGING_FEM<double>;
template void DEBUGGING_FEM<double>::Highlight_Dof<2>(BLOCK_ID,int,int,int) const;
template void DEBUGGING_FEM<double>::Highlight_Dof<3>(BLOCK_ID,int,int,int) const;
}
