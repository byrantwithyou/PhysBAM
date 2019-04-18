//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TUPLE.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include "CACHED_ELIMINATION_MATRIX.h"
#include "MATRIX_CONSTRUCTION_FEM.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_CONSTRUCTION_FEM<TV>::
MATRIX_CONSTRUCTION_FEM(COMPONENT_LAYOUT_FEM<T>& cl)
    :cl(cl)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_CONSTRUCTION_FEM<TV>::
~MATRIX_CONSTRUCTION_FEM()
{
    for(auto& p:canonical_block_matrices) delete p.key;
}

int fem_visc_table[12][12]=
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

int fem_pres_table[3][12]=
{
    {-1,0,0,1,-1,1,-1,0,0,1,1,-1},
    {0,1,0,1,-1,-1,0,0,0,2,0,-2},
    {0,0,0,2,-2,0,0,0,1,1,-1,-1}
};

//#####################################################################
// Function Fill_Canonical_Block_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Fill_Canonical_Block_Matrix(BLOCK_MATRIX<T>& mat,const CANONICAL_BLOCK<T>* cb)
{
    mat.nr=mat.nc={cb->X.m,cb->S.m,cb->X.m};
    mat.Resize();

    HASHTABLE<IV2,int> edge_lookup;
    for(int i=0;i<cb->S.m;i++)
        edge_lookup.Set(cb->S(i).Sorted(),i);

    for(IV3 v:cb->E)
    {
        IV3 dof[2]={v};
        for(int i=0;i<3;i++)
            dof[1](i)=edge_lookup.Get(v.Remove_Index(i).Sorted());
        MATRIX<T,2> F(cb->X(v.y)-cb->X(v.x),cb->X(v.z)-cb->X(v.x)),G=F.Inverse();
        T scale=mu*F.Determinant()/6;
        T p_scale=F.Determinant()/6;

        for(int r=0;r<2;r++)
            for(int s=0;s<2;s++)
                for(int i=0;i<3;i++)
                    for(int j=0;j<3;j++)
                    {
                        MATRIX<T,2> M;
                        for(int a=0;a<2;a++)
                            for(int b=0;b<2;b++)
                                M(a,b)=fem_visc_table[6*a+3*r+i][6*b+3*s+j];
                        M=scale*G.Transpose_Times(M*G);
                        mat.Add_uu(dof[r](i),r,dof[s](j),s,M+M.Trace());
                    }

        for(int s=0;s<2;s++)
            for(int i=0;i<3;i++)
                for(int j=0;j<3;j++)
                {
                    TV u;
                    for(int b=0;b<2;b++)
                        u(b)=fem_pres_table[i][6*b+3*s+j];
                    u=-p_scale*G.Transpose_Times(u);
                    mat.Add_pu(v(i),dof[s](j),s,u);
                    mat.Add_up(dof[s](j),s,v(i),u);
                }
    }
}


// entry * J / 360
int fem_u_dot_v_table[6][6]=
{
    {6,-1,-1,-4,0,0},
    {-1,6,-1,0,-4,0},
    {-1,-1,6,0,0,-4},
    {-4,0,0,32,16,16},
    {0,-4,0,16,32,16},
    {0,0,-4,16,16,32}
};

//#####################################################################
// Function Times_U_Dot_V
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Times_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& w,const BLOCK_VECTOR<T>& u) const
{
    const auto& bl=cl.blocks(b);
    const auto* cb=bl.block;
    MATRIX<T,2> M=bl.xform.M;

    HASHTABLE<IV2,int> edge_lookup;
    for(int i=0;i<cb->S.m;i++)
        edge_lookup.Set(cb->S(i).Sorted(),i);

    for(IV3 v:cb->E)
    {
        IV3 dof[2]={v};
        for(int i=0;i<3;i++)
            dof[1](i)=edge_lookup.Get(v.Remove_Index(i).Sorted());
        MATRIX<T,2> F(cb->X(v.y)-cb->X(v.x),cb->X(v.z)-cb->X(v.x));
        T scale=(M*F).Determinant()/360;

        VECTOR<TV,6> r,s;
        for(int a=0;a<2;a++) for(int i=0;i<3;i++) r(i+3*a)=u.Get_u(dof[a](i),a);
        for(int i=0;i<6;i++)
            for(int j=0;j<6;j++)
                s(i)+=r(j)*fem_u_dot_v_table[i][j];
        s*=scale;
        for(int a=0;a<2;a++) for(int i=0;i<3;i++) w.Add_u(dof[a](i),a,s(i+3*a));
    }
}

int fem_p_u_table[3][6]=
{
    {2, -1, -1, 4, 8, 8},
    {-1, 2, -1, 8, 4, 8},
    {-1, -1, 2, 8, 8, 4}
};

//#####################################################################
// Function Times_U_Dot_V
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Times_P_U(BLOCK_ID b,BLOCK_VECTOR<T>& w,const ARRAY<T>& div_v,const ARRAY<T>& div_e) const
{
    const auto& bl=cl.blocks(b);
    const auto* cb=bl.block;
    MATRIX<T,2> M=bl.xform.M;

    HASHTABLE<IV2,int> edge_lookup;
    for(int i=0;i<cb->S.m;i++)
        edge_lookup.Set(cb->S(i).Sorted(),i);

    for(IV3 v:cb->E)
    {
        IV3 dof[2]={v};
        for(int i=0;i<3;i++)
            dof[1](i)=edge_lookup.Get(v.Remove_Index(i).Sorted());
        MATRIX<T,2> F(cb->X(v.y)-cb->X(v.x),cb->X(v.z)-cb->X(v.x));
        T scale=(M*F).Determinant()/120;

        VECTOR<T,6> r;
        VECTOR<T,3> s;
        for(int i=0;i<3;i++) r(i)=div_v(dof[0](i));
        for(int i=0;i<3;i++) r(i+3)=div_e(dof[1](i));
        for(int i=0;i<3;i++)
            for(int j=0;j<6;j++)
                s(i)+=r(j)*fem_p_u_table[i][j];
        s*=scale;
        for(int i=0;i<3;i++) w.Add_p(v(i),s(i));
    }
}

int fem_line_int_u_dot_v_table[3][3] = {{4, -1, 2}, {-1, 4, 2}, {2, 2, 16}};

//#####################################################################
// Function Times_U_Dot_V
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Times_Line_Integral_U_Dot_V(BLOCK_ID b,INTERVAL<int> bc_e,BLOCK_VECTOR<T>& w,const BLOCK_VECTOR<T>& u) const
{
    const auto& bl=cl.blocks(b);
    const auto* cb=bl.block;
    MATRIX<T,2> M=bl.xform.M;
    for(int e:bc_e)
    {
        VECTOR<TV,3> r(u.Get_v(cb->S(e).x),u.Get_v(cb->S(e).y),u.Get_e(e)),s;
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                s(i)+=r(j)*fem_line_int_u_dot_v_table[i][j];

        VECTOR<TV,2> X(cb->X.Subset(cb->S(e)));
        T scale=(M*(X.x-X.y)).Magnitude()/30;
        w.Add_v(cb->S(e).x,s.x*scale);
        w.Add_v(cb->S(e).y,s.y*scale);
        w.Add_e(e,s.z*scale);
    }
}

//#####################################################################
// Function Copy_Matrix_Data
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Copy_Matrix_Data(BLOCK_MATRIX<T>& A,BLOCK_ID b,
    const DOF_PAIRS& dpa,const DOF_PAIRS& dpb,BLOCK_ID ar,BLOCK_ID ac) const
{
    const BLOCK_MATRIX<T>& B=canonical_block_matrices.Get(cl.blocks(b).block);
    MATRIX<T,2> G=cl.blocks(b).xform.M.Inverse();
    MATRIX<T,2> Ma=G*cl.blocks(ar).xform.M;
    MATRIX<T,2> Mb=G*cl.blocks(ac).xform.M;
    T sa=1/sqrt(Ma.Determinant());
    T sb=1/sqrt(Mb.Determinant());
    Ma*=sa;
    Mb*=sb;
    
    PHYSBAM_ASSERT(A.nr==dpa.num_dofs_d);
    PHYSBAM_ASSERT(B.nr==dpa.num_dofs_s);
    PHYSBAM_ASSERT(A.nc==dpb.num_dofs_d);
    PHYSBAM_ASSERT(B.nc==dpb.num_dofs_s);
    for(auto p:dpa.v)
    {
        for(auto r:dpb.v) A.Add_vv(p.x,r.x,Ma.Transpose_Times(B.Get_vv(p.y,r.y)*Mb));
        for(auto r:dpb.e) A.Add_ve(p.x,r.x,Ma.Transpose_Times(B.Get_ve(p.y,r.y)*Mb));
        for(auto r:dpb.p) A.Add_vp(p.x,r.x,Ma.Transpose_Times(B.Get_vp(p.y,r.y)*sb));
    }
    for(auto p:dpa.e)
    {
        for(auto r:dpb.v) A.Add_ev(p.x,r.x,Ma.Transpose_Times(B.Get_ev(p.y,r.y)*Mb));
        for(auto r:dpb.e) A.Add_ee(p.x,r.x,Ma.Transpose_Times(B.Get_ee(p.y,r.y)*Mb));
        for(auto r:dpb.p) A.Add_ep(p.x,r.x,Ma.Transpose_Times(B.Get_ep(p.y,r.y)*sb));
    }
    for(auto p:dpa.p)
    {
        for(auto r:dpb.v) A.Add_pv(p.x,r.x,Mb.Transpose_Times(B.Get_pv(p.y,r.y)*sa));
        for(auto r:dpb.e) A.Add_pe(p.x,r.x,Mb.Transpose_Times(B.Get_pe(p.y,r.y)*sa));
    }
}
//#####################################################################
// Function Copy_Vector_Data
//#####################################################################
// Input B should be in world space
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Copy_Vector_Data(const BLOCK_VECTOR<T>& B,BLOCK_ID b,const DOF_PAIRS& dp)
{
    BLOCK_VECTOR<T>& A=rhs_block_list(b);
    if(!A.V.m) Init_Block_Vector(A,b);

    MATRIX<T,2> M=cl.blocks(b).xform.M;
    T s=1/sqrt(M.Determinant());
    M*=s;

    PHYSBAM_ASSERT(A.n==dp.num_dofs_d);
    PHYSBAM_ASSERT(B.n==dp.num_dofs_s);
    for(auto p:dp.v) A.Add_v(p.x,M.Transpose_Times(B.Get_v(p.y)));
    for(auto p:dp.e) A.Add_e(p.x,M.Transpose_Times(B.Get_e(p.y)));
    for(auto p:dp.p) A.Add_p(p.x,s*B.Get_p(p.y));
}
//#####################################################################
// Function Fill_Block_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Fill_Block_Matrix(BLOCK_MATRIX<T>& M,const REFERENCE_BLOCK_DATA& rd)
{
    BLOCK_ID b=rd.b;
    const auto& bl=cl.blocks(b);
    Init_Block_Matrix(M,b,b);

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
            for(const auto& h:irbd.pairs)
                Copy_Matrix_Data(M,h.b,h.irreg_pairs[0],h.irreg_pairs[0],b,b);
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
Fill_Connection_Matrix(BLOCK_MATRIX<T>& M,const REFERENCE_CONNECTION_DATA& cd)
{
    Init_Block_Matrix(M,cd.b[0],cd.b[1]);
    PHYSBAM_ASSERT(cl.blocks(cd.b[0]).connections(cd.con_id[0]).master);
    auto& rd0=cl.reference_block_data(cl.blocks(cd.b[0]).ref_id);
    auto& rd1=cl.reference_block_data(cl.blocks(cd.b[1]).ref_id);
    Copy_Matrix_Data(M,cd.b[0],rd0.pairs,cd.reg_pairs[1],cd.b[0],cd.b[1]);
    Copy_Matrix_Data(M,cd.b[1],cd.reg_pairs[0],rd1.pairs,cd.b[0],cd.b[1]);
}
//#####################################################################
// Function Fill_Irregular_Connection_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Fill_Irregular_Connection_Matrix(ARRAY<BLOCK_MATRIX<T>,RID_ID>& M,const REFERENCE_IRREGULAR_DATA& ri)
{
    BLOCK_ID bb=cl.irregular_connections(ri.ic_id).regular;
    for(RID_ID j(0);j<ri.pairs.m;j++)
    {
        auto& z=ri.pairs(j);
        Init_Block_Matrix(M(j),bb,z.b);
        Copy_Matrix_Data(M(j),z.b,z.irreg_pairs[0],cl.reference_block_data(cl.blocks(z.b).ref_id).pairs,bb,z.b);
        Copy_Matrix_Data(M(j),bb,cl.reference_block_data(cl.blocks(bb).ref_id).pairs,z.irreg_pairs[1],bb,z.b);

        if(j>RID_ID(0))
        {
            RID_ID k=j-1;
            for(CON_ID cc(0);cc<cl.blocks(z.b).connections.m;cc++)
            {
                auto& c=cl.blocks(z.b).connections(cc);
                if(!c.is_regular || ri.pairs(k).b!=c.id) continue;
                const auto& dp0=cl.Regular_Connection_Pair(z.b,cc,false);
                const auto& dp1=cl.Regular_Connection_Pair(z.b,cc,true);
                Copy_Matrix_Data(M(k),z.b,z.irreg_pairs[0],dp0,bb,c.id);
                Copy_Matrix_Data(M(j),c.id,ri.pairs(k).irreg_pairs[0],dp1,bb,z.b);

                REFERENCE_BLOCK_ID ref0=cl.blocks(z.b).ref_id,ref1=cl.blocks(c.id).ref_id;
                if(c.master)
                {
                    auto key=std::make_tuple(ref0,cc,ref1,c.con_id);
                    auto& M=regular_system_blocks(cl.regular_connection_hash.Get(key));
                    Copy_Matrix_Data(M,bb,z.irreg_pairs[1],ri.pairs(k).irreg_pairs[1],z.b,c.id);
                }
                else
                {
                    auto key=std::make_tuple(ref1,c.con_id,ref0,cc);
                    auto& M=regular_system_blocks(cl.regular_connection_hash.Get(key));
                    Copy_Matrix_Data(M,bb,ri.pairs(k).irreg_pairs[1],z.irreg_pairs[1],c.id,z.b);
                }
            }
        }
    }
}
//#####################################################################
// Function Compute_Matrix_Blocks
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Compute_Matrix_Blocks()
{
    for(auto& bl:cl.blocks)
    {
        auto pr=canonical_block_matrices.Insert(bl.block,{});
        if(pr.y) Fill_Canonical_Block_Matrix(*pr.x,bl.block);
    }

    diagonal_system_blocks.Resize(cl.reference_block_data.m);
    for(REFERENCE_BLOCK_ID i(0);i<cl.reference_block_data.m;i++)
    {
        auto& rd=cl.reference_block_data(i);
        auto& M=diagonal_system_blocks(i);
        Fill_Block_Matrix(M,rd);
    }

    regular_system_blocks.Resize(cl.reference_connection_data.m);
    for(REFERENCE_CONNECTION_ID i(0);i<cl.reference_connection_data.m;i++)
    {
        auto& rc=cl.reference_connection_data(i);
        auto& M=regular_system_blocks(i);
        Fill_Connection_Matrix(M,rc);
    }

    irregular_system_blocks.Resize(cl.reference_irregular_data.m);
    for(REFERENCE_IRREGULAR_ID i(0);i<cl.reference_irregular_data.m;i++)
    {
        auto& ri=cl.reference_irregular_data(i);
        auto& is=irregular_system_blocks(i);
        is.Resize(ri.pairs.m);
        Fill_Irregular_Connection_Matrix(is,ri);
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
        CANONICAL_BLOCK<T>* cb=bl.block;
        MATRIX<T,TV::m> M=bl.xform.M.Inverse();
            
        BLOCK_VECTOR<T> w,u;
        Init_Block_Vector(w,cb);
        Init_Block_Vector(u,cb);

        for(int i=0;i<bc.bc_v.Size();i++) u.Add_v(bc.bc_v.min_corner+i,M*bc.data_v(i));
        for(int i=0;i<bc.bc_e.Size();i++) u.Add_e(bc.bc_e.min_corner+i,M*bc.data_e(i));
        w.V=canonical_block_matrices.Get(bl.block).M*-u.V;
        w.Transform(bl.xform.M,1);
        Apply_To_RHS(bc.b,w);
    }

    for(const auto& bc:cl.bc_t)
    {
        BLOCK<T>& bl=cl.blocks(bc.b);
        CANONICAL_BLOCK<T>* cb=bl.block;
        MATRIX<T,TV::m> M=bl.xform.M.Transposed()/bl.xform.M.Determinant();

        BLOCK_VECTOR<T> w,u;
        Init_Block_Vector(w,cb);
        Init_Block_Vector(u,cb);
        for(int i=0;i<bc.bc_v.Size();i++) u.Add_v(bc.bc_v.min_corner+i,M*bc.data_v(i));
        for(int i=0;i<bc.bc_e.Size();i++) u.Add_e(bc.bc_e.min_corner+i,M*bc.data_e(i));
        Times_Line_Integral_U_Dot_V(bc.b,bc.bc_e,w,u);
        w.Transform(bl.xform.M,1);
        Apply_To_RHS(bc.b,w);
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
        cem.block_list(id).M.Exchange(M.M);
    }

    regular_system_blocks.Resize(cl.reference_connection_data.m);
    for(REFERENCE_CONNECTION_ID i(0);i<cl.reference_connection_data.m;i++)
    {
        auto& M=regular_system_blocks(i);
        int id=cem.Create_Matrix_Block(false);
        reg_id(i)=id;
        cem.block_list(id).M.Exchange(M.M);
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
            cem.block_list(id).M.Exchange(is(j).M);
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
                auto i=cl.regular_connection_hash.Get(std::make_tuple(bl.ref_id,cc,cl.blocks(c.id).ref_id,c.con_id));
                int id=reg_id(i);
                cem.Add_Block_Matrix_Entry(Value(b),Value(c.id),id);
            }
        }
    }

    for(auto& ic:cl.irregular_connections)
    {
        for(RID_ID j(0);j<irreg_id(ic.ref_id).m;j++)
        {
            int id=irreg_id(ic.ref_id)(j);
            BLOCK_ID b=cl.reference_irregular_data(ic.ref_id).pairs(j).b;
            cem.Add_Block_Matrix_Entry(Value(ic.regular),Value(b),id);
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
Init_Block_Matrix(BLOCK_MATRIX<T>& M,BLOCK_ID a,BLOCK_ID b) const
{
    const auto& c=cl.reference_block_data(cl.blocks(a).ref_id);
    const auto& d=cl.reference_block_data(cl.blocks(b).ref_id);
    M.nr=c.num_dofs_d;
    M.nc=d.num_dofs_d;
    M.Resize();
}
//#####################################################################
// Function Init_Block_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Init_Block_Vector(BLOCK_VECTOR<T>& V,BLOCK_ID b) const
{
    const auto& c=cl.reference_block_data(cl.blocks(b).ref_id);
    V.n=c.num_dofs_d;
    V.Resize();
}
//#####################################################################
// Function Init_Block_Matrix
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Init_Block_Vector(BLOCK_VECTOR<T>& V,const CANONICAL_BLOCK<T>* cb) const
{
    V.n={cb->X.m,cb->S.m,cb->X.m};
    V.Resize();
}
//#####################################################################
// Function Apply_To_RHS
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Apply_To_RHS(BLOCK_ID b,const BLOCK_VECTOR<T>& w)
{
    const auto& rd=cl.reference_block_data(cl.blocks(b).ref_id);
    const auto& bl=cl.blocks(b);

    Copy_Vector_Data(w,b,rd.pairs);

    for(CON_ID cc(0);cc<bl.connections.m;cc++)
    {
        const auto& c=bl.connections(cc);
        if(c.is_regular)
            Copy_Vector_Data(w,c.id,cl.Regular_Connection_Pair(b,cc,false));
        else
        {
            const auto& ic=cl.irregular_connections(c.irreg_id);
            const auto& irbd=cl.reference_irregular_data(ic.ref_id);
            for(const auto& h:irbd.pairs)
                Copy_Vector_Data(w,h.b,h.irreg_pairs[1]);
        }
    }

    for(auto e:bl.edge_on)
    {
        const auto& ic=cl.irregular_connections(e.x);
        const auto& irbd=cl.reference_irregular_data(ic.ref_id);
        if(irbd.mapping(e.y).y)
        {
            const auto& h=irbd.pairs(irbd.mapping(e.y).x);
            Copy_Vector_Data(w,ic.regular,h.irreg_pairs[0]);
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
            Init_Block_Vector(rhs_block_list(b),b);
            continue;
        }
        auto& U=rhs_block_list(b);
        U.V=cem.vector_list(j);
        auto A = cl.blocks(b).xform.M;
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
Transform_To_World_Space(BLOCK_MATRIX<T>& M,const BLOCK_MATRIX<T>& B,BLOCK_ID a,BLOCK_ID b) const
{
    MATRIX<T,2> Ma=cl.blocks(a).xform.M.Inverse();
    MATRIX<T,2> Mb=cl.blocks(b).xform.M.Inverse();
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
Dump_Matrix_Block(SYSTEM_MATRIX_HELPER<T>& h,ARRAY<int,BLOCK_ID> first[3],const BLOCK_MATRIX<T>& M,BLOCK_ID b0,BLOCK_ID b1) const
{
    BLOCK_MATRIX<T> W;
    Init_Block_Matrix(W,b0,b1);
    Transform_To_World_Space(W,M,b0,b1);
    const auto& a=cl.reference_block_data(cl.blocks(b0).ref_id);
    const auto& b=cl.reference_block_data(cl.blocks(b1).ref_id);
    int A[3]={a.num_dofs_d.v*2,a.num_dofs_d.e*2,a.num_dofs_d.p};
    int B[3]={b.num_dofs_d.v*2,b.num_dofs_d.e*2,b.num_dofs_d.p};
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
        first[0](b)=next_u;
        next_u+=c.num_dofs_d.v*2;
        first[1](b)=next_u;
        next_u+=c.num_dofs_d.e*2;
        first[2](b)=next_p;
        next_p+=c.num_dofs_d.p;
    }
    first[2]+=next_u;
    return next_u+next_p;
}
//#####################################################################
// Function Dump_World_Space_System
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Dump_World_Space_System() const
{
    ARRAY<int,BLOCK_ID> first[3];
    int size=Compute_Global_Dof_Mapping(first);

    SYSTEM_MATRIX_HELPER<T> h;

    ARRAY<int,REFERENCE_BLOCK_ID> diag_id(cl.reference_block_data.m);
    ARRAY<int,REFERENCE_CONNECTION_ID> reg_id(cl.reference_connection_data.m);
    ARRAY<ARRAY<int,RID_ID>,REFERENCE_IRREGULAR_ID> irreg_id(cl.reference_irregular_data.m);

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        Dump_Matrix_Block(h,first,diagonal_system_blocks(cl.blocks(b).ref_id),b,b);

    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& bl=cl.blocks(b);
        for(CON_ID cc(0);cc<bl.connections.m;cc++)
        {
            const auto& c=bl.connections(cc);
            if(c.is_regular && c.master)
            {
                auto i=cl.regular_connection_hash.Get(std::make_tuple(bl.ref_id,cc,cl.blocks(c.id).ref_id,c.con_id));
                Dump_Matrix_Block(h,first,regular_system_blocks(i),b,c.id);
            }
        }
    }

    for(auto& ic:cl.irregular_connections)
    {
        for(RID_ID j(0);j<irreg_id(ic.ref_id).m;j++)
        {
            BLOCK_ID b=cl.reference_irregular_data(ic.ref_id).pairs(j).b;
            Dump_Matrix_Block(h,first,irregular_system_blocks(ic.ref_id)(j),ic.regular,b);
        }
    }

    SPARSE_MATRIX_FLAT_MXN<T> SM;
    h.Set_Matrix(size,size,SM);
    OCTAVE_OUTPUT<T>("M.txt").Write("M",SM);
}
//#####################################################################
// Function Dump_World_Space_System
//#####################################################################
template<class TV> void MATRIX_CONSTRUCTION_FEM<TV>::
Dump_World_Space_Vector(const char* name) const
{
    ARRAY<int,BLOCK_ID> first[3];
    int size=Compute_Global_Dof_Mapping(first);

    ARRAY<T> sol(size);
    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        const auto& a=cl.reference_block_data(cl.blocks(b).ref_id);
        int A[3]={a.num_dofs_d.v*2,a.num_dofs_d.e*2,a.num_dofs_d.p};
        auto& U=rhs_block_list(b);
        for(int i=0,ar=0;i<3;i++)
        {
            for(int r=0;r<A[i];r++)
                sol(first[i](b)+r)=U.V.m?U.V(r+ar):0;
            ar+=A[i];
        }
    }
    OCTAVE_OUTPUT<T>((name+(std::string)".txt").c_str()).Write(name,sol);
}
template class MATRIX_CONSTRUCTION_FEM<VECTOR<double,2> >;
}
