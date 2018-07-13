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
//#####################################################################
// Function Get_Orig_By_Blocks
//#####################################################################
template<class T> MATRIX_MXN<T>& CACHED_ELIMINATION_MATRIX<T>::
Get_Orig_By_Blocks(int b0,int b1)
{
    return *block_list(blocks_to_canonical_block_id.Get({b0,b1}));
}
template<class T,class TV>
void Setup_Block_Map(CACHED_ELIMINATION_MATRIX<T>& cem,const FLUID_LAYOUT<TV>& fl)
{
    HASHTABLE<VECTOR<int,3>,int> H;

    // matrix_id 0 is identity matrix
    // matrix_id 1 is zero matrix
    cem.block_list.Append(0);
    cem.block_list.Append(0);

    auto func=[&cem,&H,&fl](auto& a,auto& b,int s)
        {
            VECTOR<int,2> k(a.block_id,b.block_id);
            if(cem.blocks_to_canonical_block_id.Contains(k)) return;
            VECTOR<int,3> m(fl.blocks(a.block_id).block_type,fl.blocks(b.block_id).block_type,s);
            if(m.x==6 && m.y==6) assert(m.z);
            int i=-1;
            if(!H.Get(m,i)){
                i=cem.block_list.m;
                H.Set(m,cem.block_list.m);
                cem.block_list.Append(new MATRIX_MXN<T>(cem.orig_sizes(a.block_id),cem.orig_sizes(b.block_id)));}
            cem.blocks_to_canonical_block_id.Set(k,i);
        };

    for(FACE_RANGE_ITERATOR<TV::m> it(fl.used_faces.domain_indices,RF::skip_outer);it.Valid();it.Next())
    {
        auto& c0=fl.used_cells(it.face.First_Cell_Index());
        auto& c1=fl.used_cells(it.face.Second_Cell_Index());
        if(c0.global_id<0 || c1.global_id<0) continue;
        assert(c0.block_id>=0);
        assert(c1.block_id>=0);
        assert(c0.block_dof>=0);
        assert(c1.block_dof>=0);
        if(c0.block_id==c1.block_id)
            func(c0,c1,-1);
        else{
            func(c0,c1,it.face.axis*2);
            func(c1,c0,it.face.axis*2+1);}
    }
}

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
    if(rhs(r)) *rhs(r)=*block_list(inv)**rhs(r);
    for(int i=0;i<row.m;i++){
        if(row(i).c==r)
            row(i).matrix_id=id_block;
        else
            row(i).matrix_id=Compute_Mul(inv,row(i).matrix_id);}

    for(int i=0;i<row.m;i++){
        int s=row(i).c;
        if(s==r) continue;
        int elim_mat=Get_Block_Lazy(s,r);
        if(rhs(r)){
            ARRAY<T> a=*block_list(elim_mat)**rhs(r);
            if(!rhs(s)) rhs(s)=new ARRAY<T>(-a);
            else *rhs(s)-=a;}
        for(int j=0;j<row.m;j++){
            int t=row(j).c;
            int& m=Get_Block(s,t);
            m=Compute_Elim(m,elim_mat,row(j).matrix_id);}
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
// Function Compute_Inv
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Inv(int a)
{
    if(a==id_block) return id_block;
    int& r=cached_ops.Get_Or_Insert({op_inv,a,0},invalid_block);
    if(r==invalid_block){
        auto* M=new MATRIX_MXN<T>;
        block_list(a)->PLU_Inverse(*M);
        int n=block_list.Append(M);
        r=n;}
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
    int& r=cached_ops.Get_Or_Insert({op_mul,a,b},invalid_block);
    if(r==invalid_block){
        int n=block_list.Append(new MATRIX_MXN<T>(*block_list(a)**block_list(b)));
        r=n;}
    return r;
}
//#####################################################################
// Function Compute_Sub
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Sub(int a,int b)
{
    if(b==zero_block) return a;
    int& r=cached_ops.Get_Or_Insert({op_sub,a,b},invalid_block);
    if(r==invalid_block){
        auto* A=block_list(a);
        auto* B=block_list(b);
        auto* C=new MATRIX_MXN<T>(B->m,B->n);
        if(A) *C=*A-*B;
        else *C=-*B;
        int n=block_list.Append(C);
        r=n;}
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
            if(m!=zero_block) printf("%4i",m);
            else printf("    ");}
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
            if(m!=zero_block) printf("%4i",m);
            else printf("    ");}
        printf("\n");}
}
//#####################################################################
// Function Fill_Blocks
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Fill_Blocks(ARRAY<VECTOR<int,2> >& dof_map,const ARRAY<TRIPLE<int,int,T> >& data,const ARRAY<T>& rhs_vector)
{
    for(auto t:data){
        auto& M=Get_Orig_By_Blocks(dof_map(t.x).x,dof_map(t.y).x);
        M(dof_map(t.x).y,dof_map(t.y).y)=t.z;}

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
            if(e.c!=r && rhs(e.c))
                *rhs(r)-=*block_list(e.matrix_id)**rhs(e.c);}
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
                *out(r)+=*block_list(e.matrix_id)**in(e.c);}
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
template struct CACHED_ELIMINATION_MATRIX<double>;
template void Setup_Block_Map<double,VECTOR<double,2> >(
    CACHED_ELIMINATION_MATRIX<double>&,FLUID_LAYOUT<VECTOR<double,2> > const&);
}
