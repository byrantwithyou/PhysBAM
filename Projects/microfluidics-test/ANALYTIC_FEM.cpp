//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/pow.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
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
Check_Analytic_Solution(bool dump) const
{
    if(!analytic_velocity || !analytic_pressure) return false;
    T m=mc.cl.unit_m,s=mc.cl.unit_s,kg=mc.cl.unit_kg;
    T unit_p=kg*pow<2-TV::m>(m)/(s*s);
    T max_u=0,max_p=0;
    T l2_u=0,l2_p=0;
    int num_u=0,num_p=0;
    for(BLOCK_ID b(0);b<mc.cl.blocks.m;b++)
    {
        const auto& bl=mc.cl.blocks(b);
        const auto& rd=mc.cl.reference_block_data(mc.cl.blocks(b).ref_id);
        const auto& W=mc.rhs_block_list(b);
        DOF_LAYOUT<TV> dl(mc.cl,mc.cl.reference_block_data(bl.ref_id),true);
        Visit_Compressed_Dofs(dl,rd,
            [this,&bl,&W,m,s,&max_u,&l2_u,&num_u,dump](int v,const TV& Y)
            {
                TV X=xform(bl.xform,Y);
                TV U=analytic_velocity->v(X/m,0)*m/s;
                TV V=W.Get_v(v);
                max_u=std::max(max_u,(U-V).Max_Abs());
                l2_u+=(U-V).Magnitude_Squared();
                num_u++;
                if(dump)
                {
                    Add_Debug_Particle(X,VECTOR<T,3>(1,1,1));
                    Debug_Particle_Set_Attribute<TV>("V",U-V);
                }
            },
            [this,&bl,&W,m,s,&max_u,&l2_u,&num_u,dump](int e,const TV& Y)
            {
                TV X=xform(bl.xform,Y);
                TV U=analytic_velocity->v(X/m,0)*m/s;
                TV V=W.Get_e(e);
                max_u=std::max(max_u,(U-V).Max_Abs());
                l2_u+=(U-V).Magnitude_Squared();
                num_u++;
                if(dump)
                {
                    Add_Debug_Particle(X,VECTOR<T,3>(1,1,1));
                    Debug_Particle_Set_Attribute<TV>("V",U-V);
                }
            },
            [this,&bl,&W,m,unit_p,s,&max_p,&l2_p,&num_p,dump](int i,const TV& Y)
            {
                TV X=xform(bl.xform,Y);
                T p=analytic_pressure->f(X/m,0)*unit_p;
                T q=W.Get_p(i);
                max_p=std::max(max_p,std::abs(p-q));
                l2_p+=sqr(p-q);
                num_p++;
                if(dump)
                {
                    Add_Debug_Particle(X,VECTOR<T,3>(1,1,1));
                    Debug_Particle_Set_Attribute<TV>("display_size",abs(p-q));
                }
            }
        );
    }
    if(num_u) l2_u=sqrt(l2_u/num_u);
    if(num_p) l2_p=sqrt(l2_p/num_p);
    LOG::printf("u l-inf %P   u l-2 %P   p l-inf %P   p l-2 %P\n",max_u,l2_u,max_p,l2_p);
    if(dump) Flush_Frame<TV>("error");
    return max(max_u,l2_u,max_p,l2_p)<1e-10;
}
//#####################################################################
// Function Compute_RHS
//#####################################################################
template<class TV> void ANALYTIC_FEM<TV>::
Compute_RHS()
{
    T m=mc.cl.unit_m,s=mc.cl.unit_s;
    mc.rhs_block_list.Resize(mc.cl.blocks.m);

    for(BLOCK_ID b(0);b<mc.cl.blocks.m;b++)
    {
        BLOCK<T>& bl=mc.cl.blocks(b);
        CANONICAL_BLOCK<T>* cb=bl.block;
        MATRIX<T,TV::m> M=To_Dim<TV::m>(bl.xform.M.Inverse());

        BLOCK_VECTOR<TV> w,u;
        mc.Init_Block_Vector(w,b,false);
        mc.Init_Block_Vector(u,b,false);
        DOF_LAYOUT<TV> dl(mc.cl,mc.cl.reference_block_data(bl.ref_id),false);
        Visit_Wall_Dofs(dl,cb->bc_v,cb->bc_e,
            [M,&u,this,bl,m,s](const VISIT_ALL_DOFS<TV>& va)
            {
                TV Z=xform(bl.xform,va.X);
                u.Add_v(va.i,M*analytic_velocity->v(Z/m,0)*m/s);
            },
            [M,&u,this,bl,m,s](const VISIT_ALL_DOFS<TV>& va)
            {
                TV Z=xform(bl.xform,va.X);
                u.Add_e(va.i,M*analytic_velocity->v(Z/m,0)*m/s);
            });
        w.V=mc.Canonical_Matrix(b).M*-u.V;

        // TODO: handle det(M)!=1
        u.V.Fill(0);
        ARRAY<T> div_v(dl.counts.v),div_e(dl.counts.e);
        Visit_Dofs(dl,LAYER_RANGE::ALL,
            [M,&u,m,s,this,bl,&div_v](const VISIT_ALL_DOFS<TV>& va)
            {
                TV Z=xform(bl.xform,va.X);
                u.Add_v(va.i,M*Force(Z));
                div_v(va.i)=-analytic_velocity->dX(Z/m,0).Trace()/s;
            },
            [M,&u,m,s,this,bl,&div_e](const VISIT_ALL_DOFS<TV>& va)
            {
                TV Z=xform(bl.xform,va.X);
                u.Add_e(va.i,M*Force(Z));
                div_e(va.i)=-analytic_velocity->dX(Z/m,0).Trace()/s;
            });
        mc.Times_U_Dot_V(b,w,u);
        mc.Times_P_U(b,w,div_v,div_e);
        w.Transform(To_Dim<TV::m>(bl.xform.M),1);
        mc.Apply_To_RHS(b,w);
    }

    for(const auto& bc:mc.cl.bc_t)
    {
        BLOCK<T>& bl=mc.cl.blocks(bc.b);
        MATRIX<T,TV::m> M=To_Dim<TV::m>(bl.xform.M.Transposed()/bl.xform.M.Determinant());

        BLOCK_VECTOR<TV> w,u;
        mc.Init_Block_Vector(w,bc.b,false);
        mc.Init_Block_Vector(u,bc.b,false);

        DOF_LAYOUT<TV> dl(mc.cl,mc.cl.reference_block_data(bl.ref_id),false);
        Visit_Dofs<true,false>(dl,LAYER_RANGE::ALL,bc.bc_v,bc.bc_e,
            [M,&u,&bc,bl,this](const VISIT_ALL_DOFS<TV>& va)
            {
                u.Add_v(va.i,M*Traction(TV(bc.normal),xform(bl.xform,va.X)));
            },
            [M,&u,&bc,bl,this](const VISIT_ALL_DOFS<TV>& va)
            {
                u.Add_e(va.i,M*Traction(TV(bc.normal),xform(bl.xform,va.X)));
            });
        mc.Times_Line_Integral_U_Dot_V(bc.b,bc.bc_e,w,u);
        w.Transform(To_Dim<TV::m>(bl.xform.M),1);
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
    T unit_p=kg*pow<2-TV::m>(m)/(s*s);
    SYMMETRIC_MATRIX<T,TV::m> stress=analytic_velocity->dX(X/m,0).Twice_Symmetric_Part()/s*mc.mu;
    stress-=analytic_pressure->f(X/m,0)*unit_p;
    return stress*N;
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> TV ANALYTIC_FEM<TV>::
Force(const TV& X) const
{
    T m=mc.cl.unit_m,s=mc.cl.unit_s,kg=mc.cl.unit_kg;
    T unit_p=kg*pow<2-TV::m>(m)/(s*s);
    SYMMETRIC_TENSOR<T,0,TV::m> ddU=analytic_velocity->ddX(X/m,0)/(m*s);
    TV f=analytic_pressure->dX(X/m,0)*unit_p/m;
    f-=mc.mu*(Contract<1,2>(ddU)+Contract<0,2>(ddU));
    return f;
}
template struct ANALYTIC_FEM<VECTOR<double,2> >;
template struct ANALYTIC_FEM<VECTOR<double,3> >;
}
