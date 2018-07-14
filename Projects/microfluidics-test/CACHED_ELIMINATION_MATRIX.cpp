//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include "CACHED_ELIMINATION_MATRIX.h"
#include "FLUID_LAYOUT.h"

namespace PhysBAM{

template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Fill_Orig_Rows()
{
    rows.Resize(orig_sizes.m);
    for(auto h:blocks_to_canonical_block_id)
        rows(h.key.x).Append({h.key.y,h.data});

    valid_row.Resize(orig_sizes.m,use_init,true);
}

//#####################################################################
// Function Eliminate_Row
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Eliminate_Row(int r)
{
    int diag_matrix=Get_Block_Lazy(r,r);
    int inv=Compute_Inv(diag_matrix);
    ARRAY<MATRIX_BLOCK>& row=rows(r);
    if(rhs(r)) *rhs(r)=*block_list(inv).M**rhs(r);
    for(int i=0;i<row.m;i++){
        if(row(i).c==r)
            row(i).matrix_id=id_block;
        else
            row(i).matrix_id=Compute_Mul(inv,row(i).matrix_id);
        auto M=block_list(row(i).matrix_id&~use_trans).M;
        if(!M) continue;
        if(row(i).matrix_id&use_trans){PHYSBAM_ASSERT(M->Columns()==orig_sizes(r));}
        else{PHYSBAM_ASSERT(M->Rows()==orig_sizes(r));}
    }

    for(int i=0;i<row.m;i++){
        int s=row(i).c;
        if(s==r) continue;
        int elim_mat=Get_Block_Lazy(s,r);
        if(rhs(r)){
            ARRAY<T> a;
            if(elim_mat&use_trans) a=block_list(elim_mat&~use_trans).M->Transpose_Times(*rhs(r));
            else a=*block_list(elim_mat).M**rhs(r);
            if(!rhs(s)) rhs(s)=new ARRAY<T>(-a);
            else *rhs(s)-=a;}
        for(int j=0;j<row.m;j++){
            int t=row(j).c;
            int& m=Get_Block(s,t);
            m=Compute_Elim(m,elim_mat,row(j).matrix_id);
            auto M=block_list(m&~use_trans).M;
            if(!M) continue;
            if(m&use_trans){PHYSBAM_ASSERT(M->Columns()==orig_sizes(s));}
            else{PHYSBAM_ASSERT(M->Rows()==orig_sizes(s));}}
        for(int j=rows(s).m-1;j>=0;j--)
            if(rows(s)(j).matrix_id==zero_block)
                rows(s).Remove_Index_Lazy(j);}
    valid_row(r)=false;
    elimination_order.Append(r);
}
//#####################################################################
// Function Get_Block
//#####################################################################
template<class T> int& CACHED_ELIMINATION_MATRIX<T>::
Get_Block(int r,int c)
{
    ARRAY<MATRIX_BLOCK>& row=rows(r);
    for(int i=0;i<row.m;i++)
        if(row(i).c==c)
            return row(i).matrix_id;
    row.Append({c,zero_block});
    return row.Last().matrix_id;
}
//#####################################################################
// Function Get_Block
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Get_Block_Lazy(int r,int c) const
{
    const ARRAY<MATRIX_BLOCK>& row=rows(r);
    for(int i=0;i<row.m;i++)
        if(row(i).c==c)
            return row(i).matrix_id;
    return zero_block;
}
//#####################################################################
// Function Transposed
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Transposed(int a) const
{
    if(Symmetric(a)) return a;
    return a^use_trans;
}
//#####################################################################
// Function Transposed
//#####################################################################
template<class T> bool CACHED_ELIMINATION_MATRIX<T>::
Symmetric(int a) const
{
    return block_list(a&~use_trans).sym;
}
//#####################################################################
// Function Compute_Inv
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Inv(int a)
{
    if(a==id_block) return id_block;
    int& r=cached_ops.Get_Or_Insert({op_inv,a,0},invalid_block);
    if(r==invalid_block){
        auto* M=new MATRIX_MXN<T>;
        PHYSBAM_ASSERT(block_list(a).sym);
        block_list(a).M->PLU_Inverse(*M);
        int n=block_list.Append({M,true,{block_list.m}});
        r=n;

        auto ch=[this](int a){if(a&use_trans) return 't';if(Symmetric(a)) return 's';return ' ';};

        printf("inv %i%c -> %i%c\n",
            a&~use_trans,ch(a),
            r&~use_trans,ch(r)
        );
    }
    return r;
}
//#####################################################################
// Function Compute_Mul
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Mul(int a,int b)
{
    if(a==id_block) return b;
    if(b==id_block) return a;
    if(a==zero_block || b==zero_block) return zero_block;
    if((a&b&use_trans) || (!Symmetric(a) && Symmetric(b)))
        return Transposed(Compute_Mul(Transposed(b),Transposed(a)));

    int& r=cached_ops.Get_Or_Insert({op_mul,a,b},invalid_block);
    if(r==invalid_block){
        ARRAY<int> prod_list=Prod_List(a);
        prod_list.Append_Elements(Prod_List(b));
        if(int* p=prod_lookup.Get_Pointer(prod_list))
            return *p;

        auto& A=*block_list(a&~use_trans).M;
        auto& B=*block_list(b&~use_trans).M;
        MATRIX_MXN<T>& M=*new MATRIX_MXN<T>;
        if(b&use_trans) M=A.Times_Transpose(B);
        else if(a&use_trans) M=A.Transpose_Times(B);
        else M=A*B;
        ARRAY<int> prod_list_trans=Transposed(prod_list);
        bool sym=prod_list==prod_list_trans;

        int n=block_list.Append({&M,sym,prod_list});
        if(!sym) cached_ops.Set({op_mul,Transposed(b),Transposed(a)},n^use_trans);
        prod_lookup.Set(prod_list,n);
        if(!sym) prod_lookup.Set(prod_list_trans,Transposed(n));
        r=n;

        auto ch=[this](int a){if(a&use_trans) return 't';if(Symmetric(a)) return 's';return ' ';};

        printf("mul %i%c * %i%c -> %i%c\n",
            a&~use_trans,ch(a),
            b&~use_trans,ch(b),
            r&~use_trans,ch(r)
        );
    }
    return r;
}
//#####################################################################
// Function Compute_Sub
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Sub(int a,int b)
{
    if(b==zero_block) return a;
    if((a&use_trans) || (Symmetric(a) && (b&use_trans)))
        return Transposed(Compute_Sub(Transposed(a),Transposed(b)));

    int& r=cached_ops.Get_Or_Insert({op_sub,a,b},invalid_block);
    if(r==invalid_block){
        auto& A=*block_list(a&~use_trans).M;
        auto& B=*block_list(b&~use_trans).M;
        auto& C=*new MATRIX_MXN<T>;
        if(a==zero_block) C=-B;
        else if(Symmetric(a) || Symmetric(b) || !(b&use_trans)) C=A-B;
        else C=A-B.Transposed();
        int n=block_list.Append({&C,Symmetric(a)&&Symmetric(b),{block_list.m}});
        r=n;

        auto ch=[this](int a){if(a&use_trans) return 't';if(Symmetric(a)) return 's';return ' ';};
    
        printf("sub %i%c - %i%c -> %i%c\n",
            a&~use_trans,ch(a),
            b&~use_trans,ch(b),
            r&~use_trans,ch(r)
        );
    }
    return r;
}
//#####################################################################
// Function Compute_Elim
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Elim(int a,int b,int c)
{
    if(a==b && c==id_block) return zero_block;
    return Compute_Sub(a,Compute_Mul(b,c));
}
//#####################################################################
// Function Print_Full
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Print_Full() const
{
    for(int i=0;i<rows.m;i++){
        for(int j=0;j<rows.m;j++){
            int m=Get_Block_Lazy(i,j);
            if(m==zero_block) printf("    ");
            else if(m&use_trans) printf("%3it",m&~use_trans);
            else printf("%4i",m);}
        printf("\n");}
}
//#####################################################################
// Function Print_Current
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Print_Current() const
{
    for(int i=0;i<rows.m;i++){
        if(!valid_row(i)) continue;
        for(int j=0;j<rows.m;j++){
            if(!valid_row(j)) continue;
            int m=Get_Block_Lazy(i,j);
            if(m==zero_block) printf("    ");
            else if(m&use_trans) printf("%3it",m&~use_trans);
            else printf("%4i",m);}
        printf("\n");}
}
//#####################################################################
// Function Fill_Blocks
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Fill_Blocks(ARRAY<VECTOR<int,2> >& dof_map,const ARRAY<VECTOR<int,3> >& coded_entries,
    const ARRAY<T>& code_values,const ARRAY<T>& rhs_vector)
{
    HASHTABLE<VECTOR<int,2>,ARRAY<VECTOR<int,3> > > coded_blocks;
    for(auto t:coded_entries){
        ARRAY<VECTOR<int,3> >& a=coded_blocks.Get_Or_Insert({dof_map(t.x).x,dof_map(t.y).x});
        a.Append({dof_map(t.x).y,dof_map(t.y).y,t.z});}
    
    block_list.Append({0,true,{zero_block}});
    block_list.Append({0,true,{id_block}});
    HASHTABLE<int,int> hash_to_id;
    for(auto& t:coded_blocks){
        t.data.Sort(LEXICOGRAPHIC_COMPARE());
        int h=Hash(t.data),id=-1;
        if(!hash_to_id.Get(h,id)){
            id=block_list.m;
            hash_to_id.Set(h,id);
            MATRIX_MXN<T>* M=new MATRIX_MXN<T>(orig_sizes(t.key.x),orig_sizes(t.key.y));
            block_list.Append({M,t.key.x==t.key.y,{block_list.m}});
            for(auto& s:t.data){
                (*M)(s.x,s.y)=code_values(s.z);
                std::swap(s.x,s.y);}
            t.data.Sort(LEXICOGRAPHIC_COMPARE());
            int th=Hash(t.data);
            if(th!=h)
                hash_to_id.Set(th,id|use_trans);}
        blocks_to_canonical_block_id.Set(t.key,id);}

    Unpack_Vector(dof_map,rhs,rhs_vector);
}
//#####################################################################
// Function Back_Solve
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Back_Solve()
{
    for(int j=elimination_order.m-1;j>=0;j--){
        int r=elimination_order(j);
        if(!rhs(r))
            rhs(r)=new ARRAY<T>(orig_sizes(r));
        for(auto e:rows(r))
            if(e.c!=r && rhs(e.c)){
                if(e.matrix_id&use_trans) *rhs(r)-=block_list(e.matrix_id).M->Transpose_Times(*rhs(e.c));
                else *rhs(r)-=*block_list(e.matrix_id).M**rhs(e.c);}}
}
//#####################################################################
// Function Test_State
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Add_Times(ARRAY<ARRAY<T>*>& out,const ARRAY<ARRAY<T>*>& in)
{
    for(int r=0;r<rows.m;r++)
        for(auto e:rows(r))
            if(in(e.c)){
                if(!out(r)) out(r)=new ARRAY<T>(orig_sizes(r));
                if(e.matrix_id&use_trans) *out(r)+=block_list(e.matrix_id&~use_trans).M->Transpose_Times(*in(e.c));
                else *out(r)+=*block_list(e.matrix_id).M**in(e.c);}
}
//#####################################################################
// Function Test_State
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Test_State()
{
    LOG::printf("BEGIN TEST\n");
    ARRAY<ARRAY<T>*> residual(test_sol.m);
    for(int i=0;i<test_sol.m;i++){
        if(rhs(i)) residual(i)=new ARRAY<T>(-*rhs(i));
        else residual(i)=new ARRAY<T>(orig_sizes(i));}

    Add_Times(residual,test_sol);

    for(int r=0;r<rows.m;r++)
    {
        LOG::printf("A %P\n",*residual(r));
        if(rhs(r)) LOG::printf("B %P\n",*rhs(r));
    }
    LOG::printf("END TEST\n");
}
//#####################################################################
// Function Unpack_Vector
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Unpack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<ARRAY<T>*>& u,const ARRAY<T>& v)
{
    u.Delete_Pointers_And_Clean_Memory();
    u.Resize(orig_sizes.m);
    for(int i=0;i<v.m;i++)
        if(v(i)){
            int b=dof_map(i).x;
            ARRAY<T>*& a=u(b);
            if(!a) a=new ARRAY<T>(orig_sizes(b));
            (*a)(dof_map(i).y)=v(i);}
}
//#####################################################################
// Function Pack_Vector
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Pack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<T>& v,const ARRAY<ARRAY<T>*>& u)
{
    v.Resize(dof_map.m,init_all,0);
    for(int i=0;i<v.m;i++){
        int b=dof_map(i).x;
        if(u(b)) v(i)=(*u(b))(dof_map(i).y);}
}
//#####################################################################
// Function Transposed
//#####################################################################
template<class T> ARRAY<int> CACHED_ELIMINATION_MATRIX<T>::
Transposed(const ARRAY<int>& a) const
{
    ARRAY<int> r;
    for(int i=a.m-1;i>=0;i--)
        r.Append(Transposed(a(i)));
    return r;
}
//#####################################################################
// Function Prod_List
//#####################################################################
template<class T> ARRAY<int> CACHED_ELIMINATION_MATRIX<T>::
Prod_List(int a) const
{
    const ARRAY<int>& ar=block_list(a&~use_trans).prod_list;
    if(a&use_trans) return Transposed(ar);
    return ar;
}
template struct CACHED_ELIMINATION_MATRIX<double>;
}
