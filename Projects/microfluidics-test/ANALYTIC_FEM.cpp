//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "ANALYTIC_FEM.h"
#include "VISITORS_FEM.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> ANALYTIC_FEM<TV>::
ANALYTIC_FEM(MATRIX_CONSTRUCTION_FEM<TV>& mc)
    :mc(mc)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ANALYTIC_FEM<TV>::
~ANALYTIC_FEM()
{
    delete analytic_velocity;
    delete analytic_pressure;
}
//#####################################################################
// Function Check_Analytic_Solution
//#####################################################################
template<class TV> bool ANALYTIC_FEM<TV>::
Check_Analytic_Solution() const
{
    if(!analytic_velocity || !analytic_pressure) return false;
    T m=mc.cl.unit_m,s=mc.cl.unit_s,kg=mc.cl.unit_kg;
    T max_u=0,max_p=0;
    T l2_u=0,l2_p=0;
    int num_u=0,num_p=0;
    for(BLOCK_ID b(0);b<mc.cl.blocks.m;b++)
    {
        const auto& bl=mc.cl.blocks(b);
        const auto* cb=bl.block;
        auto Z=[=](int i){return bl.xform*cb->X(i);};
        const auto& rd=mc.cl.reference_block_data(mc.cl.blocks(b).ref_id);
        const auto& W=mc.rhs_block_list(b);
        for(int i=0;i<cb->X.m;i++)
            if(rd.dof_map_v(i)>=0)
            {
                TV X=Z(i);
                TV U=analytic_velocity->v(X/m,0)*m/s;
                TV V=W.Get_v(rd.dof_map_v(i));
                max_u=std::max(max_u,(U-V).Max_Abs());
                l2_u+=(U-V).Magnitude_Squared();
                num_u++;
            }

        for(int i=0;i<cb->S.m;i++)
            if(rd.dof_map_e(i)>=0)
            {
                TV X=(Z(cb->S(i).x)+Z(cb->S(i).y))/2;
                TV U=analytic_velocity->v(X/m,0)*m/s;
                TV V=W.Get_e(rd.dof_map_e(i));
                max_u=std::max(max_u,(U-V).Max_Abs());
                l2_u+=(U-V).Magnitude_Squared();
                num_u++;
            }

        for(int i=0;i<cb->X.m;i++)
            if(rd.dof_map_p(i)>=0)
            {
                TV X=Z(i);
                T p=analytic_pressure->f(X/m,0)*m/s;
                T q=W.Get_p(rd.dof_map_p(i));
                max_p=std::max(max_p,std::abs(p-q));
                l2_p+=sqr(p-q);
                num_p++;
            }
    }
    if(num_u) l2_u=sqrt(l2_u/num_u);
    if(num_p) l2_p=sqrt(l2_p/num_p);
    LOG::printf("u l-inf %P   u l-2 %P   p l-inf %P   p l-2 %P\n",max_u,l2_u,max_p,l2_p);
    return max(max_u,l2_u,max_p,l2_p)<1e-10;
}
//#####################################################################
// Function Compute_RHS
//#####################################################################
template<class TV> void ANALYTIC_FEM<TV>::
Compute_RHS()
{
    T m=mc.cl.unit_m,s=mc.cl.unit_s,kg=mc.cl.unit_kg;
    mc.rhs_block_list.Resize(mc.cl.blocks.m);

    for(BLOCK_ID b(0);b<mc.cl.blocks.m;b++)
    {
        BLOCK<T>& bl=mc.cl.blocks(b);
        CANONICAL_BLOCK<T>* cb=bl.block;
        MATRIX<T,TV::m> M=bl.xform.M.Inverse();

        BLOCK_VECTOR<T> w,u;
        mc.Init_Block_Vector(w,b,false);
        mc.Init_Block_Vector(u,b,false);
        for(auto i:cb->bc_v)
        {
            TV Z=bl.xform*cb->X(i);
            u.Add_v(i,M*analytic_velocity->v(Z/m,0)*m/s);
        }
        for(auto i:cb->bc_e)
        {
            TV Z=bl.xform*cb->X.Subset(cb->S(i)).Average();
            u.Add_e(i,M*analytic_velocity->v(Z/m,0)*m/s);
        }
        w.V=mc.canonical_block_matrices.Get(bl.block).M*-u.V;

        // TODO: handle det(M)!=1
        u.V.Fill(0);
        ARRAY<T> div_v(cb->X.m),div_e(cb->S.m);
        for(int i=0;i<cb->X.m;i++)
        {
            TV Z=bl.xform*cb->X(i);
            u.Add_v(i,M*Force(Z));
            div_v(i)=-analytic_velocity->dX(Z/m,0).Trace()/s;
        }
        for(int i=0;i<cb->S.m;i++)
        {
            TV Z=bl.xform*cb->X.Subset(cb->S(i)).Average();
            u.Add_e(i,M*Force(Z));
            div_e(i)=-analytic_velocity->dX(Z/m,0).Trace()/s;
        }
        mc.Times_U_Dot_V(b,w,u);
        mc.Times_P_U(b,w,div_v,div_e);
        w.Transform(bl.xform.M,1);
        mc.Apply_To_RHS(b,w);
    }

    for(const auto& bc:mc.cl.bc_t)
    {
        BLOCK<T>& bl=mc.cl.blocks(bc.b);
        MATRIX<T,TV::m> M=bl.xform.M.Transposed()/bl.xform.M.Determinant();

        BLOCK_VECTOR<T> w,u;
        mc.Init_Block_Vector(w,bc.b,false);
        mc.Init_Block_Vector(u,bc.b,false);

        DOF_LAYOUT<TV> dl(mc.cl,mc.cl.reference_block_data(bl.ref_id),false);
        Visit_Wall_Dofs(dl,bc.bc_v,bc.bc_e,[M,&u,&bc,&bl,this](int v,const TV& X,const VECTOR<T,TV::m-1>&)
            {
                u.Add_v(v,M*Traction(bc.normal,bl.xform*X));
            },[M,&u,&bc,&bl,this](int e,const TV& X,const VECTOR<T,TV::m-1>&)
            {
                u.Add_e(e,M*Traction(bc.normal,bl.xform*X));
            });
        mc.Times_Line_Integral_U_Dot_V(bc.b,bc.bc_e,w,u);
        w.Transform(bl.xform.M,1);
        mc.Apply_To_RHS(bc.b,w);
    }
}
//#####################################################################
// Function Traction
//#####################################################################
template<class TV> TV ANALYTIC_FEM<TV>::
Traction(const TV& N,const TV& X) const
{
    T m=mc.cl.unit_m,s=mc.cl.unit_s,kg=mc.cl.unit_kg;
    SYMMETRIC_MATRIX<T,TV::m> stress=analytic_velocity->dX(X/m,0).Twice_Symmetric_Part()/s*mc.mu;
    stress-=analytic_pressure->f(X/m,0)*kg/sqr(s);
    return stress*N;
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> TV ANALYTIC_FEM<TV>::
Force(const TV& X) const
{
    T m=mc.cl.unit_m,s=mc.cl.unit_s,kg=mc.cl.unit_kg;
    SYMMETRIC_TENSOR<T,0,TV::m> ddU=analytic_velocity->ddX(X/m,0)/(m*s);
    TV f=analytic_pressure->dX(X/m,0)*kg/(m*sqr(s));
    f-=mc.mu*(Contract<1,2>(ddU)+Contract<0,2>(ddU));
    return f;
}
template class ANALYTIC_FEM<VECTOR<double,2> >;
}
