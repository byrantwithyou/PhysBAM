//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TUPLE.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include "CACHED_ELIMINATION_MATRIX.h"
#include "MATRIX_CONSTRUCTION_FEM.h"
#include "VISITORS_FEM.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_CONSTRUCTION_FEM<TV>::
MATRIX_CONSTRUCTION_FEM(COMPONENT_LAYOUT_FEM<T>& cl,const std::string& cache_pattern,int cache_size)
    :cl(cl)
{
    if(TV::m==2)
    {
        HASHTABLE<CANONICAL_BLOCK<T>*> cbs;
        for(const auto& b:cl.blocks)
            cbs.Insert(b.block);
        canonical_matrix_cache.entries.Resize(cbs.Size());
    }
    else
        canonical_matrix_cache.entries.Resize(Value(cl.reference_block_data.m));
    canonical_matrix_cache.refill_on_pagefault=true;
    canonical_matrix_cache.Init(cache_pattern,cache_size);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_CONSTRUCTION_FEM<TV>::
~MATRIX_CONSTRUCTION_FEM()
{
    for(auto& p:canonical_block_matrices) delete p.key;
}

template<int d>
struct FEM_TABLES
{
    static const constexpr int nodes=d+1;
    static const constexpr int edges=d*(d+1)/2;
    static const constexpr int surf_edges=d*(d-1)/2;
    static const constexpr int vel=nodes+edges;
    static const constexpr int surf_vel=d+surf_edges;
    const static int visc[vel*d][vel*d],visc_den=d==2?6:30;
    const static int pres[nodes][vel*d],pres_den=d==2?6:120;
    const static int u_dot_v[vel][vel],u_dot_v_den=d==2?360:2520;
    const static int p_u[nodes][vel],p_u_den=d==2?120:360;
    const static int surf_u_dot_v[surf_vel][surf_vel],surf_u_dot_v_den=d==2?30:180;
};

template<>
const int FEM_TABLES<2>::visc[12][12]=
{
    {3,1,0,0,0,-4,3,1,0,0,0,-4},
    {1,3,0,0,0,-4,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,-1,0,4,-4,0},
    {0,0,0,8,-8,0,0,4,0,4,-4,-4},
    {0,0,0,-8,8,0,-4,0,0,-4,4,4},
    {-4,-4,0,0,0,8,0,-4,0,-4,4,4},
    {3,0,1,0,-4,0,3,0,1,0,-4,0},
    {1,0,-1,4,0,-4,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,3,0,-4,0},
    {0,0,4,4,-4,-4,0,0,0,8,0,-8},
    {0,0,-4,-4,4,4,-4,0,-4,0,8,0},
    {-4,0,0,-4,4,4,0,0,0,-8,0,8}
};

template<>
const int FEM_TABLES<3>::visc[30][30]=
{
    {3,1,0,0,-4,-1,-1,1,1,0,3,1,0,0,-4,-1,-1,1,1,0,3,1,0,0,-4,-1,-1,1,1,0},
    {1,3,0,0,-4,1,1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,-3,1,3,-1,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,1,-3,-1,3,0},
    {-4,-4,0,0,8,0,0,0,0,0,-1,-3,0,0,4,4,4,-4,-4,0,-1,-3,0,0,4,4,4,-4,-4,0},
    {-1,1,0,0,0,8,4,-8,-4,0,-4,0,0,0,4,4,0,-4,0,0,-1,1,0,0,0,8,4,-8,-4,0},
    {-1,1,0,0,0,4,8,-4,-8,0,-1,1,0,0,0,4,8,-4,-8,0,-4,0,0,0,4,0,4,0,-4,0},
    {1,-1,0,0,0,-8,-4,8,4,0,1,3,0,0,-4,-4,-4,4,4,0,0,0,0,0,0,0,0,0,0,0},
    {1,-1,0,0,0,-4,-8,4,8,0,0,0,0,0,0,0,0,0,0,0,1,3,0,0,-4,-4,-4,4,4,0},
    {0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,-4,-8,4,8,0,1,-1,0,0,0,-8,-4,8,4,0},
    {3,0,1,0,-1,-4,-1,1,0,1,3,0,1,0,-1,-4,-1,1,0,1,3,0,1,0,-1,-4,-1,1,0,1},
    {1,0,-1,0,-3,0,1,3,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,1,0,3,0,1,-4,1,-1,0,-1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,1,0,-3,-1,0,3},
    {-4,0,0,0,4,4,0,-4,0,0,-1,0,1,0,8,0,4,-8,0,-4,-1,0,1,0,8,0,4,-8,0,-4},
    {-1,0,-3,0,4,4,4,-4,0,-4,-4,0,-4,0,0,8,0,0,0,0,-1,0,-3,0,4,4,4,-4,0,-4},
    {-1,0,1,0,4,0,8,-4,0,-8,-1,0,1,0,4,0,8,-4,0,-8,-4,0,0,0,0,4,4,0,0,-4},
    {1,0,3,0,-4,-4,-4,4,0,4,1,0,-1,0,-8,0,-4,8,0,4,0,0,0,0,0,0,0,0,0,0},
    {1,0,-1,0,-4,0,-8,4,0,8,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,-8,0,-4,8,0,4},
    {0,0,0,0,0,0,0,0,0,0,1,0,-1,0,-4,0,-8,4,0,8,1,0,3,0,-4,-4,-4,4,0,4},
    {3,0,0,1,-1,-1,-4,0,1,1,3,0,0,1,-1,-1,-4,0,1,1,3,0,0,1,-1,-1,-4,0,1,1},
    {1,0,0,-1,-3,1,0,0,3,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,1,0,0,-1,1,-3,0,0,-1,3,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,3,1,1,-4,0,-1,-1},
    {-4,0,0,0,4,0,4,0,-4,0,-1,0,0,1,8,4,0,0,-8,-4,-1,0,0,1,8,4,0,0,-8,-4},
    {-1,0,0,1,4,8,0,0,-4,-8,-4,0,0,0,0,4,4,0,0,-4,-1,0,0,1,4,8,0,0,-4,-8},
    {-1,0,0,-3,4,4,4,0,-4,-4,-1,0,0,-3,4,4,4,0,-4,-4,-4,0,0,-4,0,0,8,0,0,0},
    {1,0,0,-1,-4,-8,0,0,4,8,1,0,0,-1,-8,-4,0,0,8,4,0,0,0,0,0,0,0,0,0,0},
    {1,0,0,3,-4,-4,-4,0,4,4,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,-8,-4,0,0,8,4},
    {0,0,0,0,0,0,0,0,0,0,1,0,0,3,-4,-4,-4,0,4,4,1,0,0,-1,-4,-8,0,0,4,8}
};

template<>
const int FEM_TABLES<2>::pres[3][12]=
{
    {-1,0,0,1,-1,1,-1,0,0,1,1,-1},
    {0,1,0,1,-1,-1,0,0,0,2,0,-2},
    {0,0,0,2,-2,0,0,0,1,1,-1,-1}
};

template<>
const int FEM_TABLES<3>::pres[4][30]=
{
    {-3,-1,0,0,4,-4,-4,4,4,0,-3,0,-1,0,-4,4,-4,4,0,4,-3,0,0,-1,-4,-4,4,0,4,4},
    {1,3,0,0,-4,-4,-4,4,4,0,1,0,-1,0,-8,0,-4,8,0,4,1,0,0,-1,-8,-4,0,0,8,4},
    {1,-1,0,0,0,-8,-4,8,4,0,1,0,3,0,-4,-4,-4,4,0,4,1,0,0,-1,-4,-8,0,0,4,8},
    {1,-1,0,0,0,-4,-8,4,8,0,1,0,-1,0,-4,0,-8,4,0,8,1,0,0,3,-4,-4,-4,0,4,4}
};

// entry * J / 360
template<>
const int FEM_TABLES<2>::u_dot_v[6][6]=
{
    {6,-1,-1,-4,0,0},
    {-1,6,-1,0,-4,0},
    {-1,-1,6,0,0,-4},
    {-4,0,0,32,16,16},
    {0,-4,0,16,32,16},
    {0,0,-4,16,16,32}
};

template<>
const int FEM_TABLES<3>::u_dot_v[10][10]=
{
    {6,1,1,1,-4,-4,-4,-6,-6,-6},
    {1,6,1,1,-4,-6,-6,-4,-4,-6},
    {1,1,6,1,-6,-4,-6,-4,-6,-4},
    {1,1,1,6,-6,-6,-4,-6,-4,-4},
    {-4,-4,-6,-6,32,16,16,16,16,8},
    {-4,-6,-4,-6,16,32,16,16,8,16},
    {-4,-6,-6,-4,16,16,32,8,16,16},
    {-6,-4,-4,-6,16,16,8,32,16,16},
    {-6,-4,-6,-4,16,8,16,16,32,16},
    {-6,-6,-4,-4,8,16,16,16,16,32}
};

template<>
const int FEM_TABLES<2>::p_u[3][6]=
{
    {2, -1, -1, 4, 8, 8},
    {-1, 2, -1, 8, 4, 8},
    {-1, -1, 2, 8, 8, 4}
};

template<>
const int FEM_TABLES<3>::p_u[4][10]=
{
    {0,-1,-1,-1,4,4,4,2,2,2},
    {-1,0,-1,-1,4,2,2,4,4,2},
    {-1,-1,0,-1,2,4,2,4,2,4},
    {-1,-1,-1,0,2,2,4,2,4,4}
};

template<>
const int FEM_TABLES<2>::surf_u_dot_v[3][3] =
{
    {4, -1, 2},
    {-1, 4, 2},
    {2, 2, 16}
};

template<>
const int FEM_TABLES<3>::surf_u_dot_v[6][6] =
{
    {6,-1,-1,-4,0,0},
    {-1,6,-1,0,-4,0},
    {-1,-1,6,0,0,-4},
    {-4,0,0,32,16,16},
    {0,-4,0,16,32,16},
    {0,0,-4,16,16,32}
};

//#####################################################################
// Function Fill_Canonical_Block_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Fill_Canonical_Block_Matrix(BLOCK_MATRIX<TV>& mat,const REFERENCE_BLOCK_DATA& rb)
{
    Init_Block_Matrix(mat,rb.b,rb.b,false);
    mat.Resize();
    
    DOF_LAYOUT<TV> dl(cl,rb,false);
    Visit_Elements(dl,[this,&mat](const VISIT_ELEMENT_DATA<TV>& ve)
        {
            constexpr int num_u=ve.v.m+ve.e.m;
            int dof[num_u];
            for(int i=0;i<ve.v.m;i++) dof[i]=ve.v(i);
            for(int i=0;i<ve.e.m;i++) dof[ve.v.m+i]=ve.e(i);
            MATRIX<T,TV::m> F;
            for(int i=0;i<TV::m;i++) F.Set_Column(i,ve.X(i+1)-ve.X(0));
            F/=cl.unit_m;
            MATRIX<T,TV::m> G=F.Inverse();
            T scale=mu*F.Determinant()/FEM_TABLES<TV::m>::visc_den/sqr(cl.unit_m);
            T p_scale=F.Determinant()/FEM_TABLES<TV::m>::pres_den/cl.unit_m;

            for(int i=0;i<num_u;i++)
                for(int j=0;j<num_u;j++)
                {
                    MATRIX<T,TV::m> M;
                    for(int a=0;a<TV::m;a++)
                        for(int b=0;b<TV::m;b++)
                            M(a,b)=FEM_TABLES<TV::m>::visc[num_u*a+i][num_u*b+j];
                    M=scale*G.Transpose_Times(M*G);
                    mat.Add_uu(dof[i],i>=ve.v.m,dof[j],j>=ve.v.m,M+M.Trace());
                }

            for(int i=0;i<ve.v.m;i++)
                for(int j=0;j<num_u;j++)
                {
                    TV u;
                    for(int b=0;b<TV::m;b++)
                        u(b)=FEM_TABLES<TV::m>::pres[i][num_u*b+j];
                    u=-p_scale*G.Transpose_Times(u);
                    mat.Add_pu(ve.v(i),dof[j],j>=ve.v.m,u);
                    mat.Add_up(dof[j],j>=ve.v.m,ve.v(i),u);
                }
        });
}

//#####################################################################
// Function Times_U_Dot_V
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Times_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<TV>& w,const BLOCK_VECTOR<TV>& u) const
{
    const auto& bl=cl.blocks(b);
    MATRIX<T,TV::m> M=To_Dim<TV::m>(bl.xform.M);

    DOF_LAYOUT<TV> dl(cl,cl.reference_block_data(bl.ref_id),false);
    Visit_Elements(dl,[this,M,&u,&w](const VISIT_ELEMENT_DATA<TV>& ve)
        {
            constexpr int num_u=ve.v.m+ve.e.m;
            int dof[num_u];
            for(int i=0;i<ve.v.m;i++) dof[i]=ve.v(i);
            for(int i=0;i<ve.e.m;i++) dof[ve.v.m+i]=ve.e(i);
            MATRIX<T,TV::m> F;
            for(int i=0;i<TV::m;i++) F.Set_Column(i,ve.X(i+1)-ve.X(0));
            F/=cl.unit_m;
            T scale=(M*F).Determinant()/FEM_TABLES<TV::m>::u_dot_v_den;

            VECTOR<TV,num_u> r,s;
            for(int i=0;i<num_u;i++) r(i)=u.Get_u(dof[i],i>=ve.v.m);
            for(int i=0;i<num_u;i++)
                for(int j=0;j<num_u;j++)
                    s(i)+=r(j)*FEM_TABLES<TV::m>::u_dot_v[i][j];
            s*=scale;
            for(int i=0;i<num_u;i++) w.Add_u(dof[i],i>=ve.v.m,s(i));
        });
}

//#####################################################################
// Function Times_U_Dot_V
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Times_P_U(BLOCK_ID b,BLOCK_VECTOR<TV>& w,const ARRAY<T>& div_v,const ARRAY<T>& div_e) const
{
    const auto& bl=cl.blocks(b);
    MATRIX<T,TV::m> M=To_Dim<TV::m>(bl.xform.M);

    DOF_LAYOUT<TV> dl(cl,cl.reference_block_data(bl.ref_id),false);
    Visit_Elements(dl,[this,M,&w,&div_v,&div_e](const VISIT_ELEMENT_DATA<TV>& ve)
        {
            constexpr int num_u=ve.v.m+ve.e.m;
            MATRIX<T,TV::m> F;
            for(int i=0;i<TV::m;i++) F.Set_Column(i,ve.X(i+1)-ve.X(0));
            F/=cl.unit_m;
            T scale=(M*F).Determinant()/FEM_TABLES<TV::m>::p_u_den;

            VECTOR<T,num_u> r;
            VECTOR<T,ve.v.m> s;
            for(int i=0;i<ve.v.m;i++) r(i)=div_v(ve.v(i));
            for(int i=0;i<ve.e.m;i++) r(i+ve.v.m)=div_e(ve.e(i));
            for(int i=0;i<ve.v.m;i++)
                for(int j=0;j<num_u;j++)
                    s(i)+=r(j)*FEM_TABLES<TV::m>::p_u[i][j];
            s*=scale;
            for(int i=0;i<ve.v.m;i++) w.Add_p(ve.v(i),s(i));
        });
}

//#####################################################################
// Function Times_U_Dot_V
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Times_Line_Integral_U_Dot_V(BLOCK_ID b,INTERVAL<int> bc_e,BLOCK_VECTOR<TV>& w,const BLOCK_VECTOR<TV>& u) const
{
    const auto& bl=cl.blocks(b);
    MATRIX<T,TV::m> M=To_Dim<TV::m>(bl.xform.M);

    DOF_LAYOUT<TV> dl(cl,cl.reference_block_data(bl.ref_id),false);
    Visit_Faces(dl,bc_e,[this,&u,&w,M](const VISIT_FACE_DATA<TV>& vf)
        {
            constexpr int num_u=vf.v.m+vf.e.m;
            VECTOR<TV,num_u> r,s;
            for(int i=0;i<vf.v.m;i++) r(i)=u.Get_v(vf.v(i));
            for(int i=0;i<vf.e.m;i++) r(i+vf.v.m)=u.Get_e(vf.e(i));

            for(int i=0;i<num_u;i++)
                for(int j=0;j<num_u;j++)
                    s(i)+=r(j)*FEM_TABLES<TV::m>::surf_u_dot_v[i][j];

            PHYSBAM_ASSERT(abs(abs(M.Determinant())-1)<comp_tol);
            auto f=typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE(vf.X);
            T scale=f.Size()/pow<TV::m-1>(cl.unit_m)/FEM_TABLES<TV::m>::surf_u_dot_v_den/cl.unit_m;
            for(int i=0;i<vf.v.m;i++) w.Add_v(vf.v(i),s(i)*scale);
            for(int i=0;i<vf.e.m;i++) w.Add_e(vf.e(i),s(i+vf.v.m)*scale);
        });
}

//#####################################################################
// Function Copy_Matrix_Data
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Copy_Matrix_Data(BLOCK_MATRIX<TV>& A,BLOCK_ID b,
    const DOF_PAIRS& dpa,const DOF_PAIRS& dpb,BLOCK_ID ar,BLOCK_ID ac)
{
    const BLOCK_MATRIX<TV>& B=Canonical_Matrix(b);
    MATRIX<T,2> G=cl.blocks(b).xform.M.Inverse();
    MATRIX<T,TV::m> Ma=To_Dim<TV::m>(G*cl.blocks(ar).xform.M);
    MATRIX<T,TV::m> Mb=To_Dim<TV::m>(G*cl.blocks(ac).xform.M);
    T sa=1/sqrt(Ma.Determinant());
    T sb=1/sqrt(Mb.Determinant());
    Ma*=sa;
    Mb*=sb;
    
    DOF_LAYOUT<TV> dl0a(cl,cl.reference_block_data(cl.blocks(ar).ref_id),true);
    DOF_LAYOUT<TV> dl1a(cl,cl.reference_block_data(cl.blocks(b).ref_id),false);
    DOF_LAYOUT<TV> dl0b(cl,cl.reference_block_data(cl.blocks(ac).ref_id),true);
    DOF_LAYOUT<TV> dl1b(cl,cl.reference_block_data(cl.blocks(b).ref_id),false);

    PHYSBAM_ASSERT(A.nr==dl0a.counts);
    PHYSBAM_ASSERT(B.nr==dl1a.counts);
    PHYSBAM_ASSERT(A.nc==dl0b.counts);
    PHYSBAM_ASSERT(B.nc==dl1b.counts);

    Visit_Dof_Pairs(dl0a,dl1a,dpa,
        [=,&dl0b,&dl1b,&dpb,&A,&B](int dr,int sr)
        {
            Visit_Dof_Pairs(dl0b,dl1b,dpb,
                [=,&A,&B](int dc,int sc){A.Add_vv(dr,dc,Ma.Transpose_Times(B.Get_vv(sr,sc)*Mb));},
                [=,&A,&B](int dc,int sc){A.Add_ve(dr,dc,Ma.Transpose_Times(B.Get_ve(sr,sc)*Mb));},
                [=,&A,&B](int dc,int sc){A.Add_vp(dr,dc,Ma.Transpose_Times(B.Get_vp(sr,sc)*sb));});
        },
        [=,&dl0b,&dl1b,&dpb,&A,&B](int dr,int sr)
        {
            Visit_Dof_Pairs(dl0b,dl1b,dpb,
                [=,&A,&B](int dc,int sc){A.Add_ev(dr,dc,Ma.Transpose_Times(B.Get_ev(sr,sc)*Mb));},
                [=,&A,&B](int dc,int sc){A.Add_ee(dr,dc,Ma.Transpose_Times(B.Get_ee(sr,sc)*Mb));},
                [=,&A,&B](int dc,int sc){A.Add_ep(dr,dc,Ma.Transpose_Times(B.Get_ep(sr,sc)*sb));});
        },
        [=,&dl0b,&dl1b,&dpb,&A,&B](int dr,int sr)
        {
            Visit_Dof_Pairs(dl0b,dl1b,dpb,
                [=,&A,&B](int dc,int sc){A.Add_pv(dr,dc,Mb.Transpose_Times(B.Get_pv(sr,sc)*sa));},
                [=,&A,&B](int dc,int sc){A.Add_pe(dr,dc,Mb.Transpose_Times(B.Get_pe(sr,sc)*sa));},
                [=,&A,&B](int dc,int sc){});
        });
    Release_Canonical_Matrix(b);
}
//#####################################################################
// Function Copy_Vector_Data
//#####################################################################
// Input B should be in world space
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Copy_Vector_Data(const BLOCK_VECTOR<TV>& B,BLOCK_ID a,BLOCK_ID b,const DOF_PAIRS& dp)
{
    BLOCK_VECTOR<TV>& A=rhs_block_list(b);
    if(!A.V.m) Init_Block_Vector(A,b,true);

    MATRIX<T,TV::m> M=To_Dim<TV::m>(cl.blocks(b).xform.M);
    T r=1/sqrt(M.Determinant());
    M*=r;

    DOF_LAYOUT<TV> dl0(cl,cl.reference_block_data(cl.blocks(b).ref_id),true);
    DOF_LAYOUT<TV> dl1(cl,cl.reference_block_data(cl.blocks(a).ref_id),false);
    PHYSBAM_ASSERT(A.n==dl0.counts);
    PHYSBAM_ASSERT(B.n==dl1.counts);
    Visit_Dof_Pairs(dl0,dl1,dp,
        [=,&A,&B](int d,int s){A.Add_v(d,M.Transpose_Times(B.Get_v(s)));},
        [=,&A,&B](int d,int s){A.Add_e(d,M.Transpose_Times(B.Get_e(s)));},
        [=,&A,&B](int d,int s){A.Add_p(d,r*B.Get_p(s));});
}
//#####################################################################
// Function Fill_Block_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Fill_Block_Matrix(BLOCK_MATRIX<TV>& M,const REFERENCE_BLOCK_DATA& rd)
{
    BLOCK_ID b=rd.b;
    const auto& bl=cl.blocks(b);
    Init_Block_Matrix(M,b,b,true);

    Copy_Matrix_Data(M,b,rd.pairs,rd.pairs,b,b);

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
        {
            auto& p=cl.Regular_Connection_Pair(b,cc,true);
            Copy_Matrix_Data(M,c.id,p,p,b,b);
        }
        else
        {
            const auto& ic=cl.irregular_connections(c.irreg_id);
            const auto& irbd=cl.reference_irregular_data(ic.ref_id);
            for(int i=0;i<ic.edge_on.m;i++)
            {
                const auto& p=irbd.mapping(i);
                if(p.y)
                {
                    const auto& h=irbd.pairs(p.x);
                    Copy_Matrix_Data(M,ic.edge_on(i).b,h.irreg_pairs[0],h.irreg_pairs[0],b,b);
                }
            }
        }
    }

    for(auto e:bl.edge_on)
    {
        const auto& ic=cl.irregular_connections(e.x);
        const auto& irbd=cl.reference_irregular_data(ic.ref_id);
        const auto& p=irbd.mapping(e.y);
        if(p.y)
        {
            auto& h=irbd.pairs(p.x);
            Copy_Matrix_Data(M,ic.regular,h.irreg_pairs[1],h.irreg_pairs[1],b,b);
        }
    }
}
//#####################################################################
// Function Fill_Connection_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Fill_Connection_Matrix(BLOCK_MATRIX<TV>& M,const REFERENCE_CONNECTION_DATA& cd)
{
    Init_Block_Matrix(M,cd.b[0],cd.b[1],true);
    PHYSBAM_ASSERT(cl.blocks(cd.b[0]).connections(cd.con_id[0]).master);
    auto& rd0=cl.reference_block_data(cl.blocks(cd.b[0]).ref_id);
    auto& rd1=cl.reference_block_data(cl.blocks(cd.b[1]).ref_id);
    Copy_Matrix_Data(M,cd.b[0],rd0.pairs,cd.reg_pairs[1],cd.b[0],cd.b[1]);
    Copy_Matrix_Data(M,cd.b[1],cd.reg_pairs[0],rd1.pairs,cd.b[0],cd.b[1]);

    for(int k=0;k<cl.blocks(cd.b[0]).edge_on.m;k++)
    {
        const auto& eo=cl.blocks(cd.b[0]).edge_on(k);
        const auto& ic=cl.irregular_connections(eo.x);
        const auto& ri=cl.reference_irregular_data(ic.ref_id);
        for(int j=-1;j<2;j+=2)
        {
            int n=eo.y+j;
            if(n<0 || n>=ic.edge_on.m) continue;
            if(ic.edge_on(n).b==cd.b[1])
            {
                RID_ID j0=ri.mapping(eo.y).x,j1=ri.mapping(n).x;
                PHYSBAM_ASSERT(abs(Value(j0-j1))==1);
                Copy_Matrix_Data(M,ic.regular,ri.pairs(j0).irreg_pairs[1],ri.pairs(j1).irreg_pairs[1],cd.b[0],cd.b[1]);
            }
        }
    }
}
//#####################################################################
// Function Fill_Irregular_Connection_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Fill_Irregular_Connection_Matrix(BLOCK_MATRIX<TV>& M,const REFERENCE_IRREGULAR_DATA& ri,RID_ID j)
{
    BLOCK_ID bb=cl.irregular_connections(ri.ic_id).regular;
    auto& z=ri.pairs(j);
    Init_Block_Matrix(M,bb,z.b,true);
    Copy_Matrix_Data(M,z.b,z.irreg_pairs[0],cl.reference_block_data(cl.blocks(z.b).ref_id).pairs,bb,z.b);
    Copy_Matrix_Data(M,bb,cl.reference_block_data(cl.blocks(bb).ref_id).pairs,z.irreg_pairs[1],bb,z.b);

    if(j>RID_ID(0))
    {
        RID_ID k=j-1;
        for(CON_ID cc(0);cc<cl.blocks(z.b).connections.m;cc++)
        {
            auto& c=cl.blocks(z.b).connections(cc);
            if(!c.is_regular || ri.pairs(k).b!=c.id) continue;
            const auto& dp1=cl.Regular_Connection_Pair(z.b,cc,true);
            Copy_Matrix_Data(M,c.id,ri.pairs(k).irreg_pairs[0],dp1,bb,z.b);
        }
    }
    if(j<ri.pairs.m-1)
    {
        auto& y=ri.pairs(j+1);
        for(CON_ID cc(0);cc<cl.blocks(y.b).connections.m;cc++)
        {
            auto& c=cl.blocks(y.b).connections(cc);
            if(!c.is_regular || ri.pairs(j).b!=c.id) continue;
            const auto& dp0=cl.Regular_Connection_Pair(y.b,cc,false);
            Copy_Matrix_Data(M,y.b,y.irreg_pairs[0],dp0,bb,c.id);
        }
    }
}
//#####################################################################
// Function Canonical_Matrix
//#####################################################################
template<class TV> const BLOCK_MATRIX<TV>& MATRIX_CONSTRUCTION_FEM<TV>::
Canonical_Matrix(BLOCK_ID b)
{
    const auto& bl=cl.blocks(b);
    if(TV::m==2)
    {
        int i=canonical_block_matrices.Get(bl.block);
        return canonical_matrix_cache.Use(i);
    }
    else
    {
        int i=reference_matrix(bl.ref_id);
        return canonical_matrix_cache.Use(i);
    }
}
//#####################################################################
// Function Release_Canonical_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Release_Canonical_Matrix(BLOCK_ID b)
{
    const auto& bl=cl.blocks(b);
    if(TV::m==2)
    {
        int i=canonical_block_matrices.Get(bl.block);
        canonical_matrix_cache.Release(i,false);
    }
    else
    {
        int i=reference_matrix(bl.ref_id);
        canonical_matrix_cache.Release(i,false);
    }
}
//#####################################################################
// Function Compute_Matrix_Blocks
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Compute_Matrix_Blocks()
{
    if(TV::m==3) reference_matrix.Resize(cl.reference_block_data.m);
    int next=0;
    for(REFERENCE_BLOCK_ID i(0);i<cl.reference_block_data.m;i++)
    {
        auto& rd=cl.reference_block_data(i);
        auto& bl=cl.blocks(rd.b);
        if(TV::m==2)
        {
            auto pr=canonical_block_matrices.Insert(bl.block,-1);
            if(pr.y)
            {
                canonical_matrix_cache.entries(next).fill_func=[this,i](BLOCK_MATRIX<TV>& M)
                {
                    const auto& rd=cl.reference_block_data(i);
                    Fill_Canonical_Block_Matrix(M,rd);
                };
                *pr.x=next++;
            }
        }
        else
        {
            reference_matrix(i)=Value(i);
            canonical_matrix_cache.entries(Value(i)).fill_func=[this,i](BLOCK_MATRIX<TV>& M)
            {
                const auto& rd=cl.reference_block_data(i);
                Fill_Canonical_Block_Matrix(M,rd);
            };
        }
    }

    diagonal_system_blocks.Resize(cl.reference_block_data.m);
    for(REFERENCE_BLOCK_ID i(0);i<cl.reference_block_data.m;i++)
    {
        auto& M=diagonal_system_blocks(i);
        M=[this,i](MATRIX_MXN<T>& M)
            {
                BLOCK_MATRIX<TV> bm;
                auto& rd=cl.reference_block_data(i);
                Fill_Block_Matrix(bm,rd);
                M.Exchange(bm.M);
            };
    }

    regular_system_blocks.Resize(cl.reference_connection_data.m);
    for(REFERENCE_CONNECTION_ID i(0);i<cl.reference_connection_data.m;i++)
    {
        auto& M=regular_system_blocks(i);
        M=[this,i](MATRIX_MXN<T>& M)
            {
                BLOCK_MATRIX<TV> bm;
                auto& rc=cl.reference_connection_data(i);
                Fill_Connection_Matrix(bm,rc);
                M.Exchange(bm.M);
            };
    }

    irregular_system_blocks.Resize(cl.reference_irregular_data.m);
    for(REFERENCE_IRREGULAR_ID i(0);i<cl.reference_irregular_data.m;i++)
    {
        auto& is=irregular_system_blocks(i);
        auto& ri=cl.reference_irregular_data(i);
        is.Resize(ri.pairs.m);
        for(RID_ID j(0);j<ri.pairs.m;j++)
        {
            is(j)=[this,i,j](MATRIX_MXN<T>& M)
                {
                    BLOCK_MATRIX<TV> bm;
                    auto& ri=cl.reference_irregular_data(i);
                    Fill_Irregular_Connection_Matrix(bm,ri,j);
                    M.Exchange(bm.M);
                };
        }
    }
}
//#####################################################################
// Function Compute_RHS
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Compute_RHS()
{
    rhs_block_list.Resize(cl.blocks.m);

    for(const auto& bc:cl.bc_v)
    {
        BLOCK<T>& bl=cl.blocks(bc.b);
        MATRIX<T,2> M=bl.xform.M.Inverse();
            
        BLOCK_VECTOR<TV> w,u;
        Init_Block_Vector(w,bc.b,false);
        Init_Block_Vector(u,bc.b,false);

        DOF_LAYOUT<TV> dl(cl,cl.reference_block_data(bl.ref_id),false);
        TV2 A=dl.cb->X(bc.bc_v.min_corner);
        TV2 B=dl.cb->X(bc.bc_v.max_corner-1);
        T width=(B-A).Magnitude();
        T k=bc.flow_rate*6/width;
        Visit_Dofs<false,true>(dl,LAYER_RANGE::ALL,bc.bc_v,bc.bc_e,
            [k,&u,&bc,&M](const VISIT_ALL_DOFS<TV>& va)
            {
                T z=(va.uv*((T)1-va.uv)).Product();
                u.Add_v(va.i,TV(M*-k*z*bc.normal));
            },
            [k,&u,&bc,&M](const VISIT_ALL_DOFS<TV>& va)
            {
                T z=(va.uv*((T)1-va.uv)).Product();
                u.Add_e(va.i,TV(M*-k*z*bc.normal));
            });
        w.V=Canonical_Matrix(bc.b).M*-u.V;
        Release_Canonical_Matrix(bc.b);
        w.Transform(To_Dim<TV::m>(bl.xform.M),1);
        Apply_To_RHS(bc.b,w);
    }

    for(const auto& bc:cl.bc_t)
    {
        BLOCK<T>& bl=cl.blocks(bc.b);
        MATRIX<T,2> M=bl.xform.M.Transposed()/bl.xform.M.Determinant();

        BLOCK_VECTOR<TV> w,u;
        Init_Block_Vector(w,bc.b,false);
        Init_Block_Vector(u,bc.b,false);
        DOF_LAYOUT<TV> dl(cl,cl.reference_block_data(bl.ref_id),false);
        TV tr=TV(M*bc.traction);
        Visit_Dofs<false,false>(dl,LAYER_RANGE::ALL,bc.bc_v,bc.bc_e,
            [tr,&u](const VISIT_ALL_DOFS<TV>& va)
            {
                u.Add_v(va.i,tr);
            },
            [tr,&u](const VISIT_ALL_DOFS<TV>& va)
            {
                u.Add_e(va.i,tr);
            });
        Times_Line_Integral_U_Dot_V(bc.b,bc.bc_e,w,u);
        w.Transform(To_Dim<TV::m>(bl.xform.M),1);
        Apply_To_RHS(bc.b,w);
    }

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& rb=cl.reference_block_data(cl.blocks(b).ref_id);
        DOF_LAYOUT<TV> dl(cl,rb,true);
        rhs_block_list(b).n=dl.counts;
    }
}
//#####################################################################
// Function Copy_To_CEM
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Copy_To_CEM(CACHED_ELIMINATION_MATRIX<T>& cem)
{
    cem.Begin_Fill_Blocks();
    cem.rows.Resize(Value(cl.blocks.m));
    cem.rhs.Resize(Value(cl.blocks.m));

    ARRAY<int,REFERENCE_BLOCK_ID> diag_id(cl.reference_block_data.m);
    ARRAY<int,REFERENCE_CONNECTION_ID> reg_id(cl.reference_connection_data.m);
    ARRAY<ARRAY<int,RID_ID>,REFERENCE_IRREGULAR_ID> irreg_id(cl.reference_irregular_data.m);

    for(REFERENCE_BLOCK_ID i(0);i<cl.reference_block_data.m;i++)
    {
        auto& M=diagonal_system_blocks(i);
        int id=cem.Create_Matrix_Block(true);
        diag_id(i)=id;
        cem.matrix_cache.entries(id).fill_func=M;
    }

    regular_system_blocks.Resize(cl.reference_connection_data.m);
    for(REFERENCE_CONNECTION_ID i(0);i<cl.reference_connection_data.m;i++)
    {
        auto& M=regular_system_blocks(i);
        int id=cem.Create_Matrix_Block(false);
        reg_id(i)=id;
        cem.matrix_cache.entries(id).fill_func=M;
    }

    irregular_system_blocks.Resize(cl.reference_irregular_data.m);
    for(REFERENCE_IRREGULAR_ID i(0);i<cl.reference_irregular_data.m;i++)
    {
        auto& is=irregular_system_blocks(i);
        irreg_id(i).Resize(cl.reference_irregular_data(i).pairs.m);
        for(RID_ID j(0);j<irreg_id(i).m;j++)
        {
            int id=cem.Create_Matrix_Block(false);
            irreg_id(i)(j)=id;
            cem.matrix_cache.entries(id).fill_func=is(j);
        }
    }

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        int id=diag_id(cl.blocks(b).ref_id);
        cem.Add_Block_Matrix_Entry(Value(b),Value(b),id);
    }

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& bl=cl.blocks(b);
        for(CON_ID cc(0);cc<bl.connections.m;cc++)
        {
            const auto& c=bl.connections(cc);
            if(c.is_regular && c.master)
            {
                auto i=cl.regular_connection_hash.Get(cl.Regular_Connection_Key(b,cc,c.id));
                int id=reg_id(i);
                cem.Add_Block_Matrix_Entry(Value(b),Value(c.id),id);
            }
        }
    }

    for(const auto& ic:cl.irregular_connections)
    {
        const auto& irbd=cl.reference_irregular_data(ic.ref_id);
        for(int i=0;i<ic.edge_on.m;i++)
        {
            const auto& p=irbd.mapping(i);
            if(p.y)
            {
                int id=irreg_id(ic.ref_id)(p.x);
                cem.Add_Block_Matrix_Entry(Value(ic.regular),Value(ic.edge_on(i).b),id);
            }
        }
    }

    for(BLOCK_ID i(0);i<rhs_block_list.m;i++)
    {
        int j=-1;
        if(rhs_block_list(i).V.m)
            j=cem.vector_list.Append(rhs_block_list(i).V);
        cem.rhs(Value(i))=j;
    }

    cem.End_Fill_Blocks();
    cem.valid_row.Resize(Value(cl.blocks.m),use_init,true);
}
//#####################################################################
// Function Init_Block_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Init_Block_Matrix(BLOCK_MATRIX<TV>& M,BLOCK_ID a,BLOCK_ID b,bool compressed) const
{
    const auto& c=cl.reference_block_data(cl.blocks(a).ref_id);
    const auto& d=cl.reference_block_data(cl.blocks(b).ref_id);
    DOF_LAYOUT<TV> dlc(cl,c,compressed);
    DOF_LAYOUT<TV> dld(cl,d,compressed);
    M.nr=dlc.counts;
    M.nc=dld.counts;
    M.Resize();
}
//#####################################################################
// Function Init_Block_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Init_Block_Vector(BLOCK_VECTOR<TV>& V,BLOCK_ID b,bool compressed) const
{
    const auto& c=cl.reference_block_data(cl.blocks(b).ref_id);
    DOF_LAYOUT<TV> dl(cl,c,compressed);
    V.n=dl.counts;
    V.Resize();
}
//#####################################################################
// Function Apply_To_RHS
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Apply_To_RHS(BLOCK_ID b,const BLOCK_VECTOR<TV>& w)
{
    const auto& rd=cl.reference_block_data(cl.blocks(b).ref_id);
    const auto& bl=cl.blocks(b);

    Copy_Vector_Data(w,b,b,rd.pairs);

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
            Copy_Vector_Data(w,b,c.id,cl.Regular_Connection_Pair(b,cc,false));
        else
        {
            const auto& ic=cl.irregular_connections(c.irreg_id);
            const auto& irbd=cl.reference_irregular_data(ic.ref_id);
            for(int e=0;e<ic.edge_on.m;e++)
            {
                if(irbd.mapping(e).y)
                {
                    const auto& h=irbd.pairs(irbd.mapping(e).x);
                    Copy_Vector_Data(w,b,ic.edge_on(e).b,h.irreg_pairs[1]);
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
            Copy_Vector_Data(w,b,ic.regular,h.irreg_pairs[0]);
        }
    }
}
//#####################################################################
// Function Transform_Solution
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Transform_Solution(const CACHED_ELIMINATION_MATRIX<T>& cem,bool inverse,bool transpose)
{
    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        int j=cem.rhs(Value(b));
        if(j<0)
        {
            Init_Block_Vector(rhs_block_list(b),b,true);
            continue;
        }
        auto& U=rhs_block_list(b);
        U.V=cem.vector_list(j);
        auto A = To_Dim<TV::m>(cl.blocks(b).xform.M);
        if(inverse) A=A.Inverse();
        T s=1/sqrt(A.Determinant());
        A*=s;
        if(transpose) A=A.Transposed();

        for(int i=0;i<U.n.v;i++) U.Set_v(i,A*U.Get_v(i));
        for(int i=0;i<U.n.e;i++) U.Set_e(i,A*U.Get_e(i));
        for(int i=0;i<U.n.p;i++) U.Set_p(i,s*U.Get_p(i));
    }
}
//#####################################################################
// Function Transform_To_World_Space
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Transform_To_World_Space(BLOCK_MATRIX<TV>& M,const BLOCK_MATRIX<TV>& B,BLOCK_ID a,BLOCK_ID b) const
{
    MATRIX<T,TV::m> Ma=To_Dim<TV::m>(cl.blocks(a).xform.M.Inverse());
    MATRIX<T,TV::m> Mb=To_Dim<TV::m>(cl.blocks(b).xform.M.Inverse());
    T sa=1/sqrt(Ma.Determinant());
    T sb=1/sqrt(Mb.Determinant());
    Ma*=sa;
    Mb*=sb;

    for(int p=0;p<B.nr.v;p++)
    {
        for(int r=0;r<B.nc.v;r++) M.Add_vv(p,r,Ma.Transpose_Times(B.Get_vv(p,r)*Mb));
        for(int r=0;r<B.nc.e;r++) M.Add_ve(p,r,Ma.Transpose_Times(B.Get_ve(p,r)*Mb));
        for(int r=0;r<B.nc.p;r++) M.Add_vp(p,r,Ma.Transpose_Times(B.Get_vp(p,r)*sb));
    }
    for(int p=0;p<B.nr.e;p++)
    {
        for(int r=0;r<B.nc.v;r++) M.Add_ev(p,r,Ma.Transpose_Times(B.Get_ev(p,r)*Mb));
        for(int r=0;r<B.nc.e;r++) M.Add_ee(p,r,Ma.Transpose_Times(B.Get_ee(p,r)*Mb));
        for(int r=0;r<B.nc.p;r++) M.Add_ep(p,r,Ma.Transpose_Times(B.Get_ep(p,r)*sb));
    }
    for(int p=0;p<B.nr.p;p++)
    {
        for(int r=0;r<B.nc.v;r++) M.Add_pv(p,r,Mb.Transpose_Times(B.Get_pv(p,r)*sa));
        for(int r=0;r<B.nc.e;r++) M.Add_pe(p,r,Mb.Transpose_Times(B.Get_pe(p,r)*sa));
    }
}
//#####################################################################
// Function Dump_Matrix_Block
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Dump_Matrix_Block(SYSTEM_MATRIX_HELPER<T>& h,ARRAY<int,BLOCK_ID> first[3],const BLOCK_MATRIX<TV>& M,BLOCK_ID b0,BLOCK_ID b1) const
{
    BLOCK_MATRIX<TV> W;
    Init_Block_Matrix(W,b0,b1,true);
    Transform_To_World_Space(W,M,b0,b1);
    const auto& a=cl.reference_block_data(cl.blocks(b0).ref_id);
    const auto& b=cl.reference_block_data(cl.blocks(b1).ref_id);
    DOF_LAYOUT<TV> dla(cl,a,true),dlb(cl,b,true);
    int A[3]={dla.counts.v*TV::m,dla.counts.e*TV::m,dla.counts.p};
    int B[3]={dlb.counts.v*TV::m,dlb.counts.e*TV::m,dlb.counts.p};
    for(int i=0,as=0;i<3;i++)
    {
        for(int j=0,bs=0;j<3;j++)
        {
            for(int r=0;r<A[i];r++)
                for(int s=0;s<B[j];s++)
                    if(W.M(r+as,s+bs))
                    {
                        h.data.Append({first[i](b0)+r,first[j](b1)+s,W.M(r+as,s+bs)});
                        if(b0!=b1)
                            h.data.Append({first[j](b1)+s,first[i](b0)+r,W.M(r+as,s+bs)});
                    }
            bs+=B[j];
        }
        as+=A[i];
    }
}
//#####################################################################
// Function Compute_Global_Dof_Mapping
//#####################################################################
template<class TV> int MATRIX_CONSTRUCTION_FEM<TV>::
Compute_Global_Dof_Mapping(ARRAY<int,BLOCK_ID> first[3]) const
{
    int next_u=0,next_p=0;
    for(int i=0;i<3;i++) first[i].Resize(cl.blocks.m);
    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& c=cl.reference_block_data(cl.blocks(b).ref_id);
        DOF_LAYOUT<TV> dl(cl,c,true);
        first[0](b)=next_u;
        next_u+=dl.counts.v*TV::m;
        first[1](b)=next_u;
        next_u+=dl.counts.e*TV::m;
        first[2](b)=next_p;
        next_p+=dl.counts.p;
    }
    first[2]+=next_u;
    return next_u+next_p;
}
//#####################################################################
// Function Dump_World_Space_System
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Dump_World_Space_System(ARRAY<int,BLOCK_ID> first[3],int size,SPARSE_MATRIX_FLAT_MXN<T>& SM) const
{
    SYSTEM_MATRIX_HELPER<T> h;

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        BLOCK_MATRIX<TV> bm;
        Init_Block_Matrix(bm,b,b,true);
        diagonal_system_blocks(cl.blocks(b).ref_id)(bm.M);
        Dump_Matrix_Block(h,first,bm,b,b);
    }

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& bl=cl.blocks(b);
        for(CON_ID cc(0);cc<bl.connections.m;cc++)
        {
            const auto& c=bl.connections(cc);
            if(c.is_regular && c.master)
            {
                auto i=cl.regular_connection_hash.Get(cl.Regular_Connection_Key(b,cc,c.id));
                BLOCK_MATRIX<TV> bm;
                Init_Block_Matrix(bm,b,c.id,true);
                regular_system_blocks(i)(bm.M);
                Dump_Matrix_Block(h,first,bm,b,c.id);
            }
        }
    }

    for(const auto& ic:cl.irregular_connections)
    {
        const auto& irbd=cl.reference_irregular_data(ic.ref_id);
        for(int i=0;i<ic.edge_on.m;i++)
        {
            const auto& p=irbd.mapping(i);
            if(p.y)
            {
                BLOCK_MATRIX<TV> bm;
                Init_Block_Matrix(bm,ic.regular,ic.edge_on(i).b,true);
                irregular_system_blocks(ic.ref_id)(p.x)(bm.M);
                Dump_Matrix_Block(h,first,bm,ic.regular,ic.edge_on(i).b);
            }
        }
    }

    h.Set_Matrix(size,size,SM);
}
//#####################################################################
// Function Dump_World_Space_System
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Dump_World_Space_System() const
{
    ARRAY<int,BLOCK_ID> first[3];
    int size=Compute_Global_Dof_Mapping(first);
    SPARSE_MATRIX_FLAT_MXN<T> SM;
    Dump_World_Space_System(first,size,SM);
    OCTAVE_OUTPUT<T>("M.txt").Write("M",SM);
}
//#####################################################################
// Function Dump_World_Space_Vector
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Dump_World_Space_Vector(ARRAY<int,BLOCK_ID> first[3],int size,ARRAY<T>& vec) const
{
    vec.Resize(size);
    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& a=cl.reference_block_data(cl.blocks(b).ref_id);
        DOF_LAYOUT<TV> dl(cl,a,true);
        int A[3]={dl.counts.v*TV::m,dl.counts.e*TV::m,dl.counts.p};
        auto& U=rhs_block_list(b);
        for(int i=0,ar=0;i<3;i++)
        {
            for(int r=0;r<A[i];r++)
                vec(first[i](b)+r)=U.V.m?U.V(r+ar):0;
            ar+=A[i];
        }
    }
}
//#####################################################################
// Function Dump_World_Space_Vector
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Dump_World_Space_Vector(const char* name) const
{
    ARRAY<int,BLOCK_ID> first[3];
    int size=Compute_Global_Dof_Mapping(first);
    ARRAY<T> sol;
    Dump_World_Space_Vector(first,size,sol);
    OCTAVE_OUTPUT<T>((name+(std::string)".txt").c_str()).Write(name,sol);
}
//#####################################################################
// Function Inverse_DOF_Lookup
//#####################################################################
template<class TV> std::tuple<BLOCK_ID,int,int,int> MATRIX_CONSTRUCTION_FEM<TV>::
Inverse_DOF_Lookup(int global_dof) const
{
    ARRAY<int,BLOCK_ID> first[3];
    Compute_Global_Dof_Mapping(first);

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& a=cl.reference_block_data(cl.blocks(b).ref_id);
        DOF_LAYOUT<TV> dl(cl,a,true);
        int A[3]={dl.counts.v*TV::m,dl.counts.e*TV::m,dl.counts.p};
        for(int i=0;i<3;i++)
            for(int r=0;r<A[i];r++)
                if(first[i](b)+r==global_dof)
                {
                    int dim=-7;
                    if(i!=2)
                    {
                        dim=r%TV::m;
                        r=r/TV::m;
                    }
                    return std::make_tuple(b,i,r,dim);
                }
    }
    PHYSBAM_ASSERT(false);
    return std::make_tuple(BLOCK_ID(-7),-7,-7,-7);
}
//#####################################################################
// Function Print_Statistics
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Print_Statistics() const
{
    HASHTABLE<CANONICAL_BLOCK<T>*> cbs;
    for(const auto& b:cl.blocks)
        cbs.Insert(b.block);
    LOG::printf("merge: %d\n",cl.num_merge);
    LOG::printf("blocks: %P\ncanonical-blocks: %P\nref-blocks: %P\nref-reg-con: %P\nref-irreg-con: %P\n",
        cl.blocks.m,cbs.Size(),cl.reference_block_data.m,cl.reference_connection_data.m,cl.reference_irregular_data.m);
    int p2=0,v2=0;
    int p3=0,v3=0;
    for(const auto& b:cl.blocks)
    {
        const auto& rb=cl.reference_block_data(b.ref_id);
        p2+=rb.num_dofs_d.p;
        v2+=rb.num_dofs_d.v+rb.num_dofs_d.e;
        if(TV::m==3)
        {
            DOF_LAYOUT<TV> dl(cl,rb,true);
            p3+=dl.counts.p;
            v3+=dl.counts.v+dl.counts.e;
        }
    }
    LOG::printf("dofs 2d p: %d, v: %d, total: %d\n",p2,v2*2,p2+v2*2);
    if(TV::m==3)
        LOG::printf("dofs 3d p: %d, v: %d, total: %d\n",p3,v3*3,p3+v3*3);
}
template class MATRIX_CONSTRUCTION_FEM<VECTOR<double,2> >;
template class MATRIX_CONSTRUCTION_FEM<VECTOR<double,3> >;
}
