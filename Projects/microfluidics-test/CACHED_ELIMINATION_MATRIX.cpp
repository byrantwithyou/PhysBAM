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
#include "FREQUENCY_TRACKER.h"
#include "JOB_SCHEDULER.h"
#include <cblas.h>
#include <lapacke.h>
#include <suitesparse/colamd.h>

namespace PhysBAM{

void Inverse(MATRIX_MXN<float>& A)
{
    ARRAY<int> idiv(A.m);
    int ret_f = LAPACKE_ssytrf( LAPACK_ROW_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_f);
    int ret_i = LAPACKE_ssytri( LAPACK_ROW_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_i);
    for(int r=0;r<A.m;r++)
        for(int c=0;c<r;c++)
            A(r,c)=A(c,r);
}

void Inverse(MATRIX_MXN<double>& A)
{
    ARRAY<int> idiv(A.m);
    int ret_f = LAPACKE_dsytrf( LAPACK_ROW_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_f);
    int ret_i = LAPACKE_dsytri( LAPACK_ROW_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_i);
    for(int r=0;r<A.m;r++)
        for(int c=0;c<r;c++)
            A(r,c)=A(c,r);
}

void Times_MM(MATRIX_MXN<float>& A,const MATRIX_MXN<float>& B,bool bt,const MATRIX_MXN<float>& C,bool ct)
{
    int m=bt?B.n:B.m;
    int k=bt?B.m:B.n;
    int n=ct?C.m:C.n;
    if(A.m!=m || A.n!=n) A.Resize(m,n);
    
    cblas_sgemm( CblasRowMajor, bt?CblasTrans:CblasNoTrans, ct?CblasTrans:CblasNoTrans,
        m,n,k,1,B.x.Get_Array_Pointer(),
        B.n, C.x.Get_Array_Pointer(), C.n,
        0, A.x.Get_Array_Pointer(), A.n);
}

void Times_MM(MATRIX_MXN<double>& A,const MATRIX_MXN<double>& B,bool bt,const MATRIX_MXN<double>& C,bool ct)
{
    int m=bt?B.n:B.m;
    int k=bt?B.m:B.n;
    int n=ct?C.m:C.n;
    if(A.m!=m || A.n!=n) A.Resize(m,n);
    
    cblas_dgemm( CblasRowMajor, bt?CblasTrans:CblasNoTrans, ct?CblasTrans:CblasNoTrans,
        m,n,k,1,B.x.Get_Array_Pointer(),
        B.n, C.x.Get_Array_Pointer(), C.n,
        0, A.x.Get_Array_Pointer(), A.n);
}

void Times_MV(ARRAY<float>& v,float a,const MATRIX_MXN<float>& M,bool t,const ARRAY<float>& u,float b)
{
    cblas_sgemv(CblasRowMajor,t?CblasTrans:CblasNoTrans,M.m,M.n,
        a, M.x.Get_Array_Pointer(), M.n, u.Get_Array_Pointer(), 1, b,
        v.Get_Array_Pointer(), 1);
}

void Times_MV(ARRAY<double>& v,double a,const MATRIX_MXN<double>& M,bool t,const ARRAY<double>& u,double b)
{
    cblas_dgemv(CblasRowMajor,t?CblasTrans:CblasNoTrans,M.m,M.n,
        a, M.x.Get_Array_Pointer(), M.n, u.Get_Array_Pointer(), 1, b,
        v.Get_Array_Pointer(), 1);
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
    PHYSBAM_ASSERT(valid_row(r));
    int diag_matrix=Get_Block_Lazy(r,r);
    int inv=Compute_Inv(diag_matrix);
    ARRAY<MATRIX_BLOCK>& row=rows(r);
    rhs(r)=Matrix_Times(inv,rhs(r));
    for(int i=0;i<row.m;i++){
        if(row(i).c==r)
            row(i).matrix_id=id_block;
        else
            row(i).matrix_id=Compute_Mul(inv,row(i).matrix_id);
        auto M=block_list(row(i).matrix_id&~use_trans).M;
        if(!M.m) continue;
        if(row(i).matrix_id&use_trans){PHYSBAM_ASSERT(M.Columns()==orig_sizes(r));}
        else{PHYSBAM_ASSERT(M.Rows()==orig_sizes(r));}
    }
    for(int i=0;i<row.m;i++){
        int s=row(i).c;
        if(s==r) continue;
        int elim_mat=Get_Block_Lazy(s,r);
        rhs(s)=Sub_Times(rhs(s),elim_mat,rhs(r));
        for(int j=0;j<row.m;j++){
            int t=row(j).c;
            int& m=Get_Block(s,t);
            m=Compute_Elim(m,elim_mat,row(j).matrix_id);
            auto M=block_list(m&~use_trans).M;
            if(!M.m) continue;
            if(m&use_trans){PHYSBAM_ASSERT(M.Columns()==orig_sizes(s));}
            else{PHYSBAM_ASSERT(M.Rows()==orig_sizes(s));}}
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
    if(int* r=cached_ops.Get_Pointer({op_inv,a,0})) return *r;
    PHYSBAM_ASSERT(block_list(a).sym);
    int n=block_list.Append({{},true,{block_list.m}});
    jobs.Append({op_inv,a,-1,-1,n,0});
    cached_ops.Set({op_inv,a,0},n);
    return n;
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

    if(int* r=cached_ops.Get_Pointer({op_mul,a,b})) return *r;
    ARRAY<int> prod_list=Prod_List(a);
    prod_list.Append_Elements(Prod_List(b));
    if(int* p=prod_lookup.Get_Pointer(prod_list))
        return *p;
    ARRAY<int> prod_list_trans=Transposed(prod_list);
    bool sym=prod_list==prod_list_trans;
    int n=block_list.Append({{},sym,prod_list});
    jobs.Append({op_mul,a,b,-1,n,0});

    if(!sym) cached_ops.Set({op_mul,Transposed(b),Transposed(a)},n^use_trans);
    prod_lookup.Set(prod_list,n);
    if(!sym) prod_lookup.Set(prod_list_trans,Transposed(n));
    cached_ops.Set({op_mul,a,b},n);
    return n;
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

    if(int* r=cached_ops.Get_Pointer({op_sub,a,b})) return *r;
    int n=block_list.Append({{},Symmetric(a)&&Symmetric(b),{block_list.m}});
    jobs.Append({op_sub,a,b,-1,n,0});

    cached_ops.Set({op_sub,a,b},n);
    return n;
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
    if(quiet) return;
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
    if(quiet) return;
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
    
    block_list.Append({{},true,{zero_block}});
    block_list.Append({{},true,{id_block}});
    HASHTABLE<int,int> hash_to_id;
    for(auto& t:coded_blocks){
        t.data.Sort(LEXICOGRAPHIC_COMPARE());
        int h=Hash(t.data),id=-1;
        if(!hash_to_id.Get(h,id)){
            id=block_list.m;
            hash_to_id.Set(h,id);
            block_list.Append({{orig_sizes(t.key.x),orig_sizes(t.key.y)},t.key.x==t.key.y,{block_list.m}});
            MATRIX_MXN<T>& M=block_list.Last().M;
            for(auto& s:t.data){
                M(s.x,s.y)=code_values(s.z);
                std::swap(s.x,s.y);}
            t.data.Sort(LEXICOGRAPHIC_COMPARE());
            int th=Hash(t.data);
            if(th!=h)
                hash_to_id.Set(th,id|use_trans);}
        blocks_to_canonical_block_id.Set(t.key,id);}

    Unpack_Vector(dof_map,orig_rhs,rhs_vector);

    rhs.Resize(orig_rhs.m,init_all,-1);
    for(int i=0;i<orig_rhs.m;i++)
        if(orig_rhs(i).m)
            rhs(i)=vector_list.Append(orig_rhs(i));
}
//#####################################################################
// Function Back_Solve
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Back_Solve()
{
    for(int j=elimination_order.m-1;j>=0;j--){
        int r=elimination_order(j);
        for(auto e:rows(r))
            if(e.c!=r)
                rhs(r)=Sub_Times(rhs(r),e.matrix_id,rhs(e.c));}
}
//#####################################################################
// Function Add_Times
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Matrix_Times(int m,int v)
{
    if(v<0 || m==zero_block) return -1;
    if(m==id_block) return v;

    int o=vector_list.Add_End();
    jobs.Append({op_Av,m,v,-1,o,10});
    return o;
}
//#####################################################################
// Function Add_Times
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Sub_Times(int out,int m,int v)
{
    if(v<0 || m==zero_block) return out;
    int o=vector_list.Add_End();
    if(m==id_block){
        if(out>=0) jobs.Append({op_vec_sub,out,v,-1,o,11});
        else jobs.Append({op_vec_neg,v,-1,-1,o,9});
        return o;}

    jobs.Append({op_u_sub_Av,out,m,v,o,13});
    return o;
}
//#####################################################################
// Function Unpack_Vector
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Unpack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<ARRAY<T> >& u,const ARRAY<T>& v)
{
    for(auto& i:u) i.Clean_Memory();
    u.Resize(orig_sizes.m);
    for(int i=0;i<v.m;i++)
        if(v(i)){
            int b=dof_map(i).x;
            ARRAY<T>& a=u(b);
            if(!a.m) a.Resize(orig_sizes(b));
            a(dof_map(i).y)=v(i);}
}
//#####################################################################
// Function Pack_Vector
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Pack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<T>& v,const ARRAY<ARRAY<T> >& u)
{
    v.Resize(dof_map.m,init_all,0);
    for(int i=0;i<v.m;i++){
        int b=dof_map(i).x;
        if(u(b).m) v(i)=u(b)(dof_map(i).y);}
}
//#####################################################################
// Function Pack_Vector
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Pack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<T>& v,const ARRAY<int>& u)
{
    v.Resize(dof_map.m,init_all,0);
    for(int i=0;i<v.m;i++){
        int b=dof_map(i).x;
        if(u(b)>=0) v(i)=vector_list(u(b))(dof_map(i).y);}
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
//#####################################################################
// Function Reduce_Rows_By_Frequency
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Reduce_Rows_By_Frequency(int begin,int end,int fill_limit)
{
    FREQUENCY_TRACKER ft;
    for(int i=begin;i<end;i++)
        if(rows(i).m<=fill_limit)
            ft.Add(Get_Block_Lazy(i,i),i);

    while(!ft.Empty())
    {
        std::vector<int> candidates,eliminated;
        std::set<int> modified;
        int id=ft.Most_Frequent(candidates);
        for(auto r:candidates)
            if(modified.find(r)==modified.end()){
                eliminated.push_back(r);
                ft.Remove(id,r);
                for(auto a:rows(r))
                    if(a.c!=r)
                        modified.insert(a.c);}

        for(auto r:modified)
            ft.Remove(Get_Block_Lazy(r,r),r);

        for(auto r:eliminated)
            Eliminate_Row(r);

        for(auto r:modified)
            if(r>=begin && r<end && rows(r).m<=fill_limit)
                ft.Add(Get_Block_Lazy(r,r),r);
    }
}
//#####################################################################
// Function Full_Reordered_Elimination
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Full_Reordered_Elimination()
{
    ARRAY<int> reduce(rows.m,use_init,-1),remaining_dofs;
    for(int i=0;i<reduce.m;i++)
        if(valid_row(i))
            reduce(i)=remaining_dofs.Append(i);
    if(remaining_dofs.m){
        ARRAY<int> indices,offsets;
        offsets.Append(indices.m);
        for(int i=0;i<remaining_dofs.m;i++){
            for(const auto& r:rows(remaining_dofs(i)))
                indices.Append(reduce(r.c));
            std::sort(&indices(offsets.Last()),&indices(0)+indices.m);
            offsets.Append(indices.m);}

        int stats[COLAMD_STATS];
        ARRAY<int> perm(remaining_dofs.m+1);
        int ret=symamd(remaining_dofs.m,indices.Get_Array_Pointer(),
            offsets.Get_Array_Pointer(),perm.Get_Array_Pointer(),
            0,stats,calloc,free);
        PHYSBAM_ASSERT(ret);
        perm.Pop();
        for(int i=0;i<perm.m;i++)
            Eliminate_Row(remaining_dofs(perm(i)));}
}
//#####################################################################
// Function Execute_Jobs
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Execute_Jobs(int num_threads)
{
    JOB_SCHEDULER<JOB,CACHED_ELIMINATION_MATRIX<T> > scheduler(this);
    ARRAY<int> provides_mat(block_list.m,use_init,-1);
    ARRAY<int> provides_vec(vector_list.m,use_init,-1);
    for(auto& j:jobs){
        j.job_id=scheduler.Add_Job(&j,j.priority);
        j.priority=0;
        if(j.o>=0){
            if(j.is_vec_mask&8) provides_vec(j.o)=j.job_id;
            else provides_mat(j.o)=j.job_id;}}
    for(auto j:jobs){
        int x[3]={j.a&~use_trans,j.b&~use_trans,j.c&~use_trans};
        for(int i=0;i<3;i++)
            if(x[i]>=0){
                int k=(j.is_vec_mask&(1<<i))?provides_vec(x[i]):provides_mat(x[i]);
                PHYSBAM_ASSERT(k<jobs.m);
                if(k>=0)
                    scheduler.Register_Dependency(k,j.job_id);}}
    
    scheduler.Compute_Priority_By_Paths();
    scheduler.Execute_Jobs(num_threads);
}
//#####################################################################
// Function Execute
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::JOB::
Execute(CACHED_ELIMINATION_MATRIX<T>* cem)
{
    auto& bl=cem->block_list;
    auto& vl=cem->vector_list;
    switch(op)
    {
        case op_inv:
            bl(o).M=bl(a).M;
            Inverse(bl(o).M);
            break;
        case op_mul:
            Times_MM(bl(o).M,bl(a&~use_trans).M,a&use_trans,bl(b&~use_trans).M,b&use_trans);
            break;
        case op_sub:
            {
                auto& A=bl(a&~use_trans).M;
                auto& B=bl(b&~use_trans).M;
                auto& C=bl(o).M;
                if(a==zero_block) C=-B;
                else if(cem->Symmetric(a) || cem->Symmetric(b) || !(b&use_trans)) C=A-B;
                else C=A-B.Transposed();
            }
            break;
        case op_Av:
            {
                auto& M=bl(a&~use_trans).M;
                vl(o).Resize(a&use_trans?M.n:M.m);
                Times_MV(vl(o),1,M,a&use_trans,vl(b),0);
            }
            break;
        case op_vec_neg:
            vl(o)=-vl(a);
            break;
        case op_vec_sub:
            vl(o)=vl(a)-vl(b);
            break;
        case op_u_sub_Av:
            {
                auto& M=bl(b&~use_trans).M;
                if(a<0) vl(o).Resize(b&use_trans?M.n:M.m);
                else vl(o)=vl(a);
                Times_MV(vl(o),-1,M,b&use_trans,vl(c),1);
            }
            break;
        case free_mat:
            bl(a).M.x.Clean_Memory();
            break;
        case free_vec:
            vl(a).Clean_Memory();
            break;
        default: PHYSBAM_FATAL_ERROR();
    }
}

template struct CACHED_ELIMINATION_MATRIX<double>;
}
