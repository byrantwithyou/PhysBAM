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
namespace
{
// const char* op_names[]={
//     "op_nop",
//     "op_mat_inv","op_mat_mul","op_mat_add",
//     "op_vec_mul","op_vec_add",
//     "free_mat","free_vec"};

enum op_type
{
    op_nop,
    op_mat_inv,op_mat_mul,op_mat_add,
    op_vec_mul,op_vec_add,
    free_mat,free_vec,
    op_last
};

enum {zero_block=0,id_block=1,invalid_block=-1,
      use_trans=1<<30,use_neg=1<<29,is_vec=1<<28,
      raw_mask=~(use_trans|use_neg)};

int arg_type[op_last][4] =
{
    [op_nop]={0,0,0,0},
    [op_mat_inv]={1,0,0,1},
    [op_mat_mul]={1,1,1,1},
    [op_mat_add]={1,1,0,1},
    [op_vec_mul]={2,1,2,2},
    [op_vec_add]={2,2,0,2},
    [free_mat]={1,0,0,0},
    [free_vec]={2,0,0,0}
};
}

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

// A = sa*A+sbc*op(B,bt)*op(C,ct)
void Times_MM(MATRIX_MXN<float>& A,float sa,const MATRIX_MXN<float>& B,bool bt,const MATRIX_MXN<float>& C,bool ct,float sbc)
{
    int m=bt?B.n:B.m;
    int k=bt?B.m:B.n;
    int n=ct?C.m:C.n;
    if(A.m!=m || A.n!=n) A.Resize(m,n);
    
    cblas_sgemm( CblasRowMajor, bt?CblasTrans:CblasNoTrans, ct?CblasTrans:CblasNoTrans,
        m,n,k,sbc,B.x.Get_Array_Pointer(),
        B.n, C.x.Get_Array_Pointer(), C.n,
        sa, A.x.Get_Array_Pointer(), A.n);
}

void Times_MM(MATRIX_MXN<double>& A,double sa,const MATRIX_MXN<double>& B,bool bt,const MATRIX_MXN<double>& C,bool ct,double sbc)
{
    int m=bt?B.n:B.m;
    int k=bt?B.m:B.n;
    int n=ct?C.m:C.n;
    if(A.m!=m || A.n!=n) A.Resize(m,n);
    
    cblas_dgemm( CblasRowMajor, bt?CblasTrans:CblasNoTrans, ct?CblasTrans:CblasNoTrans,
        m,n,k,sbc,B.x.Get_Array_Pointer(),
        B.n, C.x.Get_Array_Pointer(), C.n,
        sa, A.x.Get_Array_Pointer(), A.n);
}

// v = b*v + a*M*u
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
            row(i).matrix_id=Compute_Mul(inv,row(i).matrix_id);}
    for(int i=0;i<row.m;i++){
        int s=row(i).c;
        if(s==r) continue;
        int elim_mat=Get_Block_Lazy(s,r);
        rhs(s)=Sub_Times(rhs(s),elim_mat,rhs(r));
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
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Negate(int a) const
{
    if(a==zero_block) return a;
    return a^use_neg;
}
//#####################################################################
// Function Transposed
//#####################################################################
template<class T> bool CACHED_ELIMINATION_MATRIX<T>::
Symmetric(int a) const
{
    return block_list(a&raw_mask).sym;
}
//#####################################################################
// Function Compute_Inv
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Inv(int a)
{
    if(a==id_block) return id_block;
    if(a&use_neg) return Negate(Compute_Inv(Negate(a)));
    if(int* r=cached_ops.Get_Pointer({op_mat_inv,a,0})) return *r;
    PHYSBAM_ASSERT(block_list(a).sym);
    int n=block_list.Append({{},true,{block_list.m}});
    jobs.Append({op_mat_inv,{a},n});
    cached_ops.Set({op_mat_inv,a,0},n);
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
    if(a&use_neg) return Negate(Compute_Mul(Negate(a),b));
    if(b&use_neg) return Negate(Compute_Mul(a,Negate(b)));

    if(int* r=cached_ops.Get_Pointer({op_mat_mul,a,b})) return *r;
    ARRAY<int> prod_list=Prod_List(a);
    prod_list.Append_Elements(Prod_List(b));
    if(int* p=prod_lookup.Get_Pointer(prod_list))
        return *p;
    ARRAY<int> prod_list_trans=Transposed(prod_list);
    bool sym=prod_list==prod_list_trans;
    int n=block_list.Append({{},sym,prod_list});
    jobs.Append({op_mat_mul,{-1,a,b},n,{0,1}});

    if(!sym) cached_ops.Set({op_mat_mul,Transposed(b),Transposed(a)},n^use_trans);
    prod_lookup.Set(prod_list,n);
    if(!sym) prod_lookup.Set(prod_list_trans,Transposed(n));
    cached_ops.Set({op_mat_mul,a,b},n);
    return n;
}
//#####################################################################
// Function Compute_Add
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Add(int a,int b)
{
    if(b==zero_block) return a;
    if(a==zero_block) return b;
    if(a==Negate(b)) return zero_block;
    int ra=a&raw_mask,rb=b&raw_mask;
    if(ra>rb) return Compute_Add(b,a);
    if((a&use_trans) || (Symmetric(a) && (b&use_trans)))
        return Transposed(Compute_Add(Transposed(a),Transposed(b)));
    if(a&use_neg) return Negate(Compute_Add(Negate(a),Negate(b)));

    if(int* r=cached_ops.Get_Pointer({op_mat_add,a,b})) return *r;
    int n=block_list.Append({{},Symmetric(a)&&Symmetric(b),{block_list.m}});

    jobs.Append({op_mat_add,{a,b},n});
    cached_ops.Set({op_mat_add,a,b},n);
    return n;
}
//#####################################################################
// Function Compute_Elim
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Elim(int a,int b,int c)
{
    if(a==b && c==id_block) return zero_block;
    return Compute_Add(a,Negate(Compute_Mul(b,c)));
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
    jobs.Append({op_vec_mul,{-1,m,v},o,{0,1}});
    return o;
}
//#####################################################################
// Function Add_Times
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Sub_Times(int out,int m,int v)
{
    if(v<0 || m==zero_block) return out;
    if(m==id_block && out<0) return Negate(v);
    int o=vector_list.Add_End();
    if(m==id_block) jobs.Append({op_vec_add,{out,v},o});
    else jobs.Append({op_vec_mul,{out,m,v},o,{1,-1}});
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
        if(u(b)>=0){
            v(i)=vector_list(u(b)&raw_mask)(dof_map(i).y);
            if(u(b)&use_neg) v(i)=-v(i);}}
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
    Compute_Job_Deps();
    Eliminate_Trans();
    Combine_Ops();
    Simplify_Jobs();

    JOB_SCHEDULER<JOB,CACHED_ELIMINATION_MATRIX<T> > scheduler(this);
    for(int i=0;i<jobs.m;i++){
        auto& j=jobs(i);
        PHYSBAM_ASSERT(i==scheduler.Add_Job(&j,0));}
    for(int i=0;i<jobs.m;i++){
        auto& j=jobs(i);
        for(int l=0;l<3;l++)
            if(j.dep_jobs[l]>=0)
                scheduler.Register_Dependency(j.dep_jobs[l],i);}
    
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
        case op_nop: break;
        case op_mat_inv:
            bl(o).M=bl(a[0]).M;
            Inverse(bl(o).M);
            break;
        case op_mat_mul:
            {
                auto& A=bl(a[1]&raw_mask).M;
                auto& B=bl(a[2]&raw_mask).M;
                auto& C=bl(o).M;
                if(a[0]>=0) C=bl(a[0]).M;
                else C.Resize(a[1]&use_trans?A.n:A.m,a[2]&use_trans?B.m:B.n);
                Times_MM(C,s[0],A,a[1]&use_trans,B,a[2]&use_trans,s[1]);
            }
            break;
        case op_mat_add:
            {
                auto& A=bl(a[0]&raw_mask).M;
                auto& B=bl(a[1]&raw_mask).M;
                auto& C=bl(o).M;
                if(a[1]&use_trans)
                {
                    if(a[1]&use_neg) C=A-B.Transposed();
                    else C=A+B.Transposed();
                }
                else
                {
                    if(a[1]&use_neg) C=A-B;
                    else C=A+B;
                }
            }
            break;
        case op_vec_mul:
            {
                auto& M=bl(a[1]&raw_mask).M;
                if(a[0]>=0) vl(o)=vl(a[0]);
                else vl(o).Resize(a[1]&use_trans?M.n:M.m);
                Times_MV(vl(o),s[1],M,a[1]&use_trans,vl(a[2]),s[0]);
            }
            break;
        case op_vec_add:
            if(a[1]&use_neg) vl(o)=vl(a[0])-vl(a[1]);
            else vl(o)=vl(a[0])+vl(a[1]);
            break;
        case free_mat:
            bl(a[0]).M.x.Clean_Memory();
            break;
        case free_vec:
            vl(a[0]).Clean_Memory();
            break;
        default: PHYSBAM_FATAL_ERROR();
    }
}
//#####################################################################
// Function Compute_Job_Deps
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Compute_Job_Deps()
{
    ARRAY<int> provides_mat(block_list.m,use_init,-1);
    ARRAY<int> provides_vec(vector_list.m,use_init,-1);
    for(int i=0;i<jobs.m;i++){
        auto& j=jobs(i);
        if(j.o>=0){
            if(arg_type[j.op][3]==2) provides_vec(j.o)=i;
            else if(arg_type[j.op][3]==1) provides_mat(j.o)=i;}}
    for(int i=0;i<jobs.m;i++){
        auto& j=jobs(i);
        for(int l=0;l<3;l++){
            j.dep_jobs[l]=-1;
            if(arg_type[j.op][l]==0) continue;
            if(j.a[l]<0) continue;
            int k=-1;
            int a=j.a[l]&raw_mask;
            if(arg_type[j.op][l]==1) k=provides_mat(a);
            else k=provides_vec(a);
            if(k<0) continue;
            PHYSBAM_ASSERT(k<jobs.m);
            jobs(k).users.Append(i);
            j.dep_jobs[l]=k;}}
}
//#####################################################################
// Function Eliminate_Trans
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Eliminate_Trans()
{
    for(int i=0;i<jobs.m;i++)
    {
        auto& j=jobs(i);
        if(j.op==op_mat_add && (j.a[1]&use_trans))
        {
            if(j.dep_jobs[1]>=0)
            {
                auto& k=jobs(j.dep_jobs[1]);
                if(k.users.m==1)
                {
                    if(k.op==op_mat_mul)
                    {
                        std::swap(k.a[1],k.a[2]);
                        std::swap(k.dep_jobs[1],k.dep_jobs[2]);
                        k.a[1]=Transposed(k.a[1]);
                        k.a[2]=Transposed(k.a[2]);
                        j.a[1]=Transposed(j.a[1]);
                        continue;
                    }
                }
            }
        }
    }
}
//#####################################################################
// Function Compute_Job_Deps
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Simplify_Jobs()
{
    for(auto& j:jobs)
    {
        switch(j.op)
        {
            case op_mat_mul:
            case op_vec_mul:
                if(j.a[0]>=0 && (j.a[0]&use_neg))
                {
                    j.a[0]&=~use_neg;
                    j.s[0]=-j.s[0];
                }
                if(j.a[1]&use_neg)
                {
                    j.a[1]&=~use_neg;
                    j.s[1]=-j.s[1];
                }
                if(j.a[2]&use_neg)
                {
                    j.a[2]&=~use_neg;
                    j.s[1]=-j.s[1];
                }
                break;
        }
    }
}
//#####################################################################
// Function Combine_Ops
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Combine_Ops()
{
    for(int i=0;i<jobs.m;i++)
    {
        auto& j=jobs(i);
        if(j.op==op_mat_mul && j.a[0]<0 && j.users.m==1)
        {
            auto& k=jobs(j.users(0));
            if(k.op==op_mat_add)
            {
                if(k.a[1]&use_trans)
                {
                    printf("trans discard %i -> %i\n",i,j.users(0));
                    continue;
                }
                if(k.a[0]==j.o)
                {
                    j.s[0]=(k.a[1]&use_neg)?-1:1;
                    j.a[0]=k.a[1]&raw_mask;
                    j.dep_jobs[0]=k.dep_jobs[1];
                }
                else
                {
                    if(k.a[1]&use_neg) j.s[1]=-j.s[1];
                    j.s[0]=1;
                    j.a[0]=k.a[0];
                    j.dep_jobs[0]=k.dep_jobs[0];
                }
                Combine_With_Next_Job(i);
            }
        }
    }
}
//#####################################################################
// Function Combine_With_Next_Job
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Combine_With_Next_Job(int j)
{
    auto& k=jobs(j);
    int l=k.users(0);
    auto& m=jobs(l);
    k.o=m.o;
    k.users=m.users;
    for(auto& u:m.users)
    {
        for(auto&r:jobs(u).dep_jobs)
            if(r==l)
                r=j;
    }
    m.op=op_nop;
    for(auto& i:m.a) i=-1;
    for(auto& i:m.dep_jobs) i=-1;
    m.o=-1;
    m.users.Clean_Memory();
}
template struct CACHED_ELIMINATION_MATRIX<double>;
}
