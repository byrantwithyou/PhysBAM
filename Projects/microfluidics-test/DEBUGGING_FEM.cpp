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
#include "DEBUGGING_FEM.h"
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
        VECTOR<TV,3> P;
        for(int i=0;i<3;i++) P(i)=bl.xform*cb->X(t(i));
        for(auto p:P) Add_Debug_Object(VECTOR<TV,2>(p,P.Average()),VECTOR<T,3>(.5,.5,.5));
    }
    HASHTABLE<int> he,hv;
    he.Set_All(cb->bc_e);
    hv.Set_All(cb->bc_v);
    for(int i=0;i<cb->S.m;i++)
    {
        TV X=Z(cb->S(i).x);
        TV Y=Z(cb->S(i).y);
        Add_Debug_Object(VECTOR<TV,2>(X,Y),he.Contains(i)?VECTOR<T,3>(0,1,1):VECTOR<T,3>(0,0,1));
    }
    for(int i=0;i<cb->X.m;i++)
    {
        Add_Debug_Particle(Z(i),hv.Contains(i)?VECTOR<T,3>(1,1,0):VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>("display_size",.1);
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
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                },
                [&](int a,int b,bool o)
                {
                    Add_Debug_Particle((Z(cb->S(a).x)+Z(cb->S(a).y))/2,o?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                });
        }
        else
        {
            const auto& ic=cl.irregular_connections(c.irreg_id);
            Visit_Irregular_Cross_Section_Dofs(cl.blocks,ic,
                [&](const IRREGULAR_VISITOR& iv)
                {
                    TV A=Z(iv.r0),B=Z(iv.r1);
                    Add_Debug_Particle(A,VECTOR<T,3>(iv.b0,0,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                    Add_Debug_Particle(B,VECTOR<T,3>(iv.b1,0,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.25);
                    Add_Debug_Particle((A+B)/2,VECTOR<T,3>(iv.be,0,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                });
        }
    }

    for(auto i:bl.edge_on)
    {
        auto& ic=cl.irregular_connections(i.x);
        auto& id=ic.edge_on(i.y);
        bool o0=i.y>=(ic.edge_on.m/2+1);
        bool o1=(i.y>=ic.edge_on.m/2);

        TV A=Z(id.v0),B=Z(id.v1);
        Add_Debug_Particle(A,VECTOR<T,3>(1,1,1)/(1+o0));
        Debug_Particle_Set_Attribute<TV>("display_size",.05);
        Add_Debug_Particle(B,VECTOR<T,3>(1,1,1)/(1+o1));
        Debug_Particle_Set_Attribute<TV>("display_size",.06);
        Add_Debug_Particle((A+B)/2,VECTOR<T,3>(1,1,1)/(1+o1));
        Debug_Particle_Set_Attribute<TV>("display_size",.05);
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
        VECTOR<TV,3> P;
        for(int i=0;i<3;i++) P(i)=bl.xform*cb->X(t(i));
        for(auto p:P) Add_Debug_Object(VECTOR<TV,2>(p,P.Average()),VECTOR<T,3>(.5,.5,.5));
    }
    HASHTABLE<int> he,hv;
    he.Set_All(cb->bc_e);
    hv.Set_All(cb->bc_v);
    for(int i=0;i<cb->S.m;i++)
    {
        TV X=Z(cb->S(i).x);
        TV Y=Z(cb->S(i).y);
        Add_Debug_Object(VECTOR<TV,2>(X,Y),he.Contains(i)?VECTOR<T,3>(0,1,1):VECTOR<T,3>(0,0,1));
    }
    Flush_Frame<TV>(LOG::sprintf("block %P\n",b).c_str());
    for(int i=0;i<cb->X.m;i++)
    {
        Add_Debug_Text(Z(i),LOG::sprintf("%d",i),VECTOR<T,3>(1,1,0));
    }
    for(int i=0;i<cb->S.m;i++)
    {
        Add_Debug_Text((Z(cb->S(i).x)+Z(cb->S(i).y))/2,LOG::sprintf("%d",i),VECTOR<T,3>(0,1,0));
    }
    Flush_Frame<TV>(LOG::sprintf("block %P full\n",b).c_str());
    for(int i=0;i<cb->X.m;i++)
    {
        int k=rd.dof_map_v(i);
        if(k>=0) Add_Debug_Text(Z(i),LOG::sprintf("%d",k),VECTOR<T,3>(1,1,0));
    }
    Flush_Frame<TV>(LOG::sprintf("block %P v\n",b).c_str());
    for(int i=0;i<cb->S.m;i++)
    {
        int k=rd.dof_map_e(i);
        if(k>=0) Add_Debug_Text((Z(cb->S(i).x)+Z(cb->S(i).y))/2,LOG::sprintf("%d",k),VECTOR<T,3>(0,1,0));
    }
    Flush_Frame<TV>(LOG::sprintf("block %P e\n",b).c_str());
    for(int i=0;i<cb->X.m;i++)
    {
        int k=rd.dof_map_p(i);
        if(k>=0) Add_Debug_Text(Z(i),LOG::sprintf("%d",k),VECTOR<T,3>(1,1,0));
    }
    Flush_Frame<TV>(LOG::sprintf("block %P p\n",b).c_str());

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
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                },
                [&](int a,int b,bool o)
                {
                    Add_Debug_Particle((Z(cb->S(a).x)+Z(cb->S(a).y))/2,o?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                });
        }
        else
        {
            const auto& ic=cl.irregular_connections(c.irreg_id);
            Visit_Irregular_Cross_Section_Dofs(cl.blocks,ic,
                [&](const IRREGULAR_VISITOR& iv)
                {
                    TV A=Z(iv.r0),B=Z(iv.r1);
                    Add_Debug_Particle(A,VECTOR<T,3>(iv.b0,0,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                    Add_Debug_Particle(B,VECTOR<T,3>(iv.b1,0,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.25);
                    Add_Debug_Particle((A+B)/2,VECTOR<T,3>(iv.be,0,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",.2);
                });
        }
    }

    for(auto i:bl.edge_on)
    {
        auto& ic=cl.irregular_connections(i.x);
        auto& id=ic.edge_on(i.y);
        bool o0=i.y>=(ic.edge_on.m/2+1);
        bool o1=(i.y>=ic.edge_on.m/2);

        TV A=Z(id.v0),B=Z(id.v1);
        Add_Debug_Particle(A,VECTOR<T,3>(1,1,1)/(1+o0));
        Debug_Particle_Set_Attribute<TV>("display_size",.05);
        Add_Debug_Particle(B,VECTOR<T,3>(1,1,1)/(1+o1));
        Debug_Particle_Set_Attribute<TV>("display_size",.06);
        Add_Debug_Particle((A+B)/2,VECTOR<T,3>(1,1,1)/(1+o1));
        Debug_Particle_Set_Attribute<TV>("display_size",.05);
    }
}
//#####################################################################
// Function Visualize_Solution
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Solution(const BLOCK_VECTOR<T>& U,BLOCK_ID b,bool remap_dofs) const
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
            Debug_Particle_Set_Attribute<TV>("V",U.Get_v(k));
        }
    }

    for(int i=0;i<cb->S.m;i++)
    {
        int k=remap_dofs?rd.dof_map_e(i):i;
        if(k>=0)
        {
            Add_Debug_Particle((Z(cb->S(i).x)+Z(cb->S(i).y))/2,VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV>("V",U.Get_e(k));
        }
    }

    for(int i=0;i<cb->X.m;i++)
    {
        int k=remap_dofs?rd.dof_map_p(i):i;
        if(k>=0)
        {
            Add_Debug_Particle(Z(i),VECTOR<T,3>(0,1,0));
            Debug_Particle_Set_Attribute<TV>("display_size",U.Get_p(k));
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
            TV X=Z(t(i)),X0=Z(t((i+2)%3)),X1=Z(t((i+1)%3));
            TV x0=0.85*X+0.13*X0+0.02*X1,x1=0.85*X+0.02*X0+0.13*X1;
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
    Flush_Frame<TV>("flat velocity dofs");
    
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
    Flush_Frame<TV>("flat pressure dofs");
    

}
//#####################################################################
// Function Visualize_Tetrahedron
//#####################################################################
template<class T> void DEBUGGING_FEM<T>::
Visualize_Tetrahedron(BLOCK_ID b) const
{
    typedef VECTOR<T,3> TV;

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto* cb=cl.blocks(b).block;
        auto Z=[=](int i){return cl.blocks(b).xform*cb->X(i);};
        for(auto t:cb->E)
        {
            VECTOR<TV,3> P;
            for(int i=0;i<3;i++) P(i)=Z(t(i)).Append(-0.01);
            for(int i=0;i<3;i++) Add_Debug_Object(VECTOR<TV,2>(P(i),P((i+1)%3)),VECTOR<T,3>(.2,.2,.2));
        }
    }

    const auto& rb=cl.reference_block_data(cl.blocks(b).ref_id);
    auto Z=[=](const TV& X)
    {
        const MATRIX<T,2>& m=cl.blocks(b).xform.M;
        MATRIX<T,3> M(m.Column(0).Append(0),m.Column(1).Append(0),TV(0,0,1));
        return M*X+cl.blocks(b).xform.b.Append(0);
    };
    DOF_LAYOUT<TV> dl(cl,rb,false);
    Visit_Elements(dl,[Z](const VISIT_ELEMENT_DATA<TV>& vt)
    {
        VECTOR<TV,4> P;
        for(int i=0;i<4;i++) P(i)=Z(vt.X(i));
        Add_Debug_Object(VECTOR<TV,2>(P(0),P(1)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV,2>(P(0),P(2)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV,2>(P(0),P(3)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV,2>(P(1),P(2)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV,2>(P(2),P(3)),VECTOR<T,3>(1,1,1));
        Add_Debug_Object(VECTOR<TV,2>(P(3),P(1)),VECTOR<T,3>(1,1,1));
    });
}
template class DEBUGGING_FEM<double>;
}
