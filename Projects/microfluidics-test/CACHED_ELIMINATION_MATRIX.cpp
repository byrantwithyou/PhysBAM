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
#include "EXECUTE_HELPER.h"
#if USE_LAPACK
#include <cblas.h>
#include <lapacke.h>
#elif USE_MKL
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#endif
#include <suitesparse/colamd.h>
namespace PhysBAM{


// const char* op_names[]={
//     "op_nop",
//     "op_mat_inv","op_mat_mul","op_mat_add",
//     "op_vec_mul","op_vec_add"};

enum op_type
{
    op_nop,
    op_mat_inv,op_mat_mul,op_mat_add,
    op_vec_mul,op_vec_add,
    op_last
};

enum {pseudo_inv=1<<30,raw_op_mask=~pseudo_inv};

int arg_type[op_last][4] =
{
    [op_nop]={0,0,0,0},
    [op_mat_inv]={1,0,0,1},
    [op_mat_mul]={1,1,1,1},
    [op_mat_add]={1,1,0,1},
    [op_vec_mul]={2,1,2,2},
    [op_vec_add]={2,2,0,2}
};


void Inverse(MATRIX_MXN<float>& A)
{
    ARRAY<int> idiv(A.m);
    int ret_f = LAPACKE_ssytrf( LAPACK_COL_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_f);
    int ret_i = LAPACKE_ssytri( LAPACK_COL_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_i);
    for(int r=0;r<A.m;r++)
        for(int c=0;c<r;c++)
            A(r,c)=A(c,r);
}

void Inverse(MATRIX_MXN<double>& A)
{
    ARRAY<int> idiv(A.m);
    int ret_f = LAPACKE_dsytrf( LAPACK_COL_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_f);
    int ret_i = LAPACKE_dsytri( LAPACK_COL_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_i);
    for(int r=0;r<A.m;r++)
        for(int c=0;c<r;c++)
            A(r,c)=A(c,r);
}

void Pseudo_Inverse(MATRIX_MXN<float>& A)
{
    ARRAY<float> s(A.m),u(A.m*A.m),vt(A.m*A.m),superb(A.m);
    int ret=LAPACKE_sgesvd(LAPACK_COL_MAJOR,'S','S',A.m,A.n,A.x.Get_Array_Pointer(),A.m,
        s.Get_Array_Pointer(),u.Get_Array_Pointer(),A.m,vt.Get_Array_Pointer(),A.m,superb.Get_Array_Pointer());
    PHYSBAM_ASSERT(!ret);
    float tol=1e-8*A.m*s.Max();
    for(int i=0;i<A.m;i++)
    {
       float ss=s(i);
       if(s(i)>tol) ss=1.0/s(i);
       cblas_sscal(A.m,ss,&vt(i*A.m),1);
    }
    cblas_sgemm(CblasColMajor,CblasTrans,CblasTrans,
        A.m,A.m,A.m,1.0,vt.Get_Array_Pointer(),A.m,u.Get_Array_Pointer(),A.m,0.0,A.x.Get_Array_Pointer(),A.m);
}

static void Pseudo_Inverse(MATRIX_MXN<double>& A)
{
    ARRAY<double> s(A.m),u(A.m*A.m),vt(A.m*A.m),superb(A.m);
    int ret=LAPACKE_dgesvd(LAPACK_COL_MAJOR,'S','S',A.m,A.n,A.x.Get_Array_Pointer(),A.m,
        s.Get_Array_Pointer(),u.Get_Array_Pointer(),A.m,vt.Get_Array_Pointer(),A.m,superb.Get_Array_Pointer());
    PHYSBAM_ASSERT(!ret);
    double tol=1e-16*A.m*s.Max();
    for(int i=0;i<A.m;i++)
    {
       double ss=s(i);
       if(s(i)>tol) ss=1.0/s(i);
       cblas_dscal(A.m,ss,&vt(i*A.m),1);
    }
    cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,
        A.m,A.m,A.m,1.0,vt.Get_Array_Pointer(),A.m,u.Get_Array_Pointer(),A.m,0.0,A.x.Get_Array_Pointer(),A.m);
}

// A = sa*A+sbc*op(B,bt)*op(C,ct)
void Times_MM(MATRIX_MXN<float>& A,float sa,const MATRIX_MXN<float>& B,bool bt,const MATRIX_MXN<float>& C,bool ct,float sbc)
{
    int m=bt?B.n:B.m;
    int k=bt?B.m:B.n;
    int n=ct?C.m:C.n;
    if(A.m!=m || A.n!=n) A.Resize(m,n);
    
    cblas_sgemm( CblasColMajor, bt?CblasTrans:CblasNoTrans, ct?CblasTrans:CblasNoTrans,
        m,n,k,sbc,B.x.Get_Array_Pointer(),
        B.m, C.x.Get_Array_Pointer(), C.m,
        sa, A.x.Get_Array_Pointer(), A.m);
}

void Times_MM(MATRIX_MXN<double>& A,double sa,const MATRIX_MXN<double>& B,bool bt,const MATRIX_MXN<double>& C,bool ct,double sbc)
{
    int m=bt?B.n:B.m;
    int k=bt?B.m:B.n;
    int n=ct?C.m:C.n;
    if(A.m!=m || A.n!=n) A.Resize(m,n);
    
    cblas_dgemm( CblasColMajor, bt?CblasTrans:CblasNoTrans, ct?CblasTrans:CblasNoTrans,
        m,n,k,sbc,B.x.Get_Array_Pointer(),
        B.m, C.x.Get_Array_Pointer(), C.m,
        sa, A.x.Get_Array_Pointer(), A.m);
}

// v = b*v + a*M*u
void Times_MV(ARRAY<float>& v,float a,const MATRIX_MXN<float>& M,bool t,const ARRAY<float>& u,float b)
{
    cblas_sgemv(CblasColMajor,t?CblasTrans:CblasNoTrans,M.m,M.n,
        a, M.x.Get_Array_Pointer(), M.m, u.Get_Array_Pointer(), 1, b,
        v.Get_Array_Pointer(), 1);
}

void Times_MV(ARRAY<double>& v,double a,const MATRIX_MXN<double>& M,bool t,const ARRAY<double>& u,double b)
{
    cblas_dgemv(CblasColMajor,t?CblasTrans:CblasNoTrans,M.m,M.n,
        a, M.x.Get_Array_Pointer(), M.m, u.Get_Array_Pointer(), 1, b,
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
    int inv=Compute_Inv(diag_matrix,pinv_last_blk && elimination_order.m==valid_row.m-1);
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
// Function Compute_Inv
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Compute_Inv(int a,bool pinv)
{
    if(a==id_block) return id_block;
    if(a&use_neg) return Negate(Compute_Inv(Negate(a),pinv));
    int op_inv=op_mat_inv;
    if(pinv) op_inv|=pseudo_inv;
    if(int* r=cached_ops.Get_Pointer({op_inv,a,0})) return *r;
    PHYSBAM_ASSERT(block_list(a).sym);
    int n=Create_Matrix_Block(true);
    jobs.Append({op_inv,{a},n});
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
    if(a&use_neg) return Negate(Compute_Mul(Negate(a),b));
    if(b&use_neg) return Negate(Compute_Mul(a,Negate(b)));

    if(int* r=cached_ops.Get_Pointer({op_mat_mul,a,b})) return *r;
    ARRAY<int> prod_list=Prod_List(a);
    prod_list.Append_Elements(Prod_List(b));
    if(int* p=prod_lookup.Get_Pointer(prod_list))
        return *p;
    ARRAY<int> prod_list_trans=Transposed(prod_list);
    bool sym=prod_list==prod_list_trans;
    matrix_cache.entries.Add_End();
    int n=block_list.Append({sym,prod_list});
    jobs.Append({op_mat_mul,{-1,a,b},n});

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
    int n=Create_Matrix_Block(Symmetric(a)&&Symmetric(b));

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
            bool tr=m&use_trans;
            m&=~use_trans;
            if(m&use_neg) m=-(m&~use_neg);
            if(m==zero_block) printf("    ");
            else if(tr) printf("%3it",m);
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
            bool tr=m&use_trans;
            m&=~use_trans;
            if(m&use_neg) m=-(m&~use_neg);
            if(m==zero_block) printf("    ");
            else if(tr) printf("%3it",m);
            else printf("%4i",m);}
        printf("\n");}
}
//#####################################################################
// Function Begin_Fill_Blocks
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Begin_Fill_Blocks()
{
    matrix_cache.entries.Add_End();
    matrix_cache.entries.Add_End();
    block_list.Append({true,{zero_block}});
    block_list.Append({true,{id_block}});
}
//#####################################################################
// Function End_Fill_Blocks
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
End_Fill_Blocks()
{
    num_orig_blocks=block_list.m;
    num_orig_vectors=vector_list.m;
}
//#####################################################################
// Function Fill_Blocks
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Fill_Blocks(ARRAY<VECTOR<int,2>,DOF_ID>& dof_map,const ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    const ARRAY<T,CODE_ID>& code_values,const ARRAY<T,DOF_ID>& rhs_vector)
{
    Begin_Fill_Blocks();
    HASHTABLE<VECTOR<int,2>,ARRAY<TRIPLE<int,int,CODE_ID> > > coded_blocks;
    for(auto t:coded_entries){
        auto& a=coded_blocks.Get_Or_Insert({dof_map(t.x).x,dof_map(t.y).x});
        a.Append({dof_map(t.x).y,dof_map(t.y).y,t.z});}
    
    HASHTABLE<int,int> hash_to_id;
    for(auto& t:coded_blocks){
        t.data.Sort();
        int h=Hash(t.data),id=-1;
        if(!hash_to_id.Get(h,id)){
            id=block_list.m;
            hash_to_id.Set(h,id);
            int n=Create_Matrix_Block(t.key.x==t.key.y);
            matrix_cache.entries(n).fill_func=[this,&t,&code_values](MATRIX_MXN<T>& M)
                {
                    M.Resize(orig_sizes(t.key.x),orig_sizes(t.key.y));
                    for(auto& s:t.data) M(s.x,s.y)=code_values(s.z);
                };
            auto S=t.data;
            for(auto& s:S) std::swap(s.x,s.y);
            S.Sort();
            int th=Hash(S);
            if(th!=h)
                hash_to_id.Set(th,id|use_trans);
        }
        blocks_to_canonical_block_id.Set(t.key,id);}

    Unpack_Vector(dof_map,orig_rhs,rhs_vector);

    rhs.Resize(orig_rhs.m,init_all,-1);
    for(int i=0;i<orig_rhs.m;i++)
        if(orig_rhs(i).m)
            rhs(i)=vector_list.Append(orig_rhs(i));
    End_Fill_Blocks();
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
    jobs.Append({op_vec_mul,{-1,m,v},o});
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
    else jobs.Append({op_vec_mul,{out,m,Negate(v)},o});
    return o;
}
//#####################################################################
// Function Unpack_Vector
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Unpack_Vector(ARRAY<VECTOR<int,2>,DOF_ID>& dof_map,ARRAY<ARRAY<T> >& u,const ARRAY<T,DOF_ID>& v)
{
    for(auto& i:u) i.Clean_Memory();
    u.Resize(orig_sizes.m);
    for(DOF_ID i(0);i<v.m;i++)
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
Pack_Vector(ARRAY<VECTOR<int,2>,DOF_ID>& dof_map,ARRAY<T,DOF_ID>& v,const ARRAY<ARRAY<T> >& u)
{
    v.Resize(dof_map.m,init_all,0);
    for(DOF_ID i(0);i<v.m;i++){
        int b=dof_map(i).x;
        if(u(b).m) v(i)=u(b)(dof_map(i).y);}
}
//#####################################################################
// Function Pack_Vector
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Pack_Vector(ARRAY<VECTOR<int,2>,DOF_ID>& dof_map,ARRAY<T,DOF_ID>& v,const ARRAY<int>& u)
{
    v.Resize(dof_map.m,init_all,0);
    for(DOF_ID i(0);i<v.m;i++){
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
    Relabel();

    JOB_SCHEDULER<JOB,CACHED_ELIMINATION_MATRIX<T> > scheduler(this);
    for(int i=0;i<jobs.m;i++){
        int id=scheduler.Add_Job(&jobs(i),0);
        PHYSBAM_ASSERT(i==id);}
    for(auto a:dep_list)
        scheduler.Register_Dependency(a.x,a.y);

    scheduler.Compute_Priority_By_Paths();
    scheduler.Execute_Jobs(num_threads);
}
//#####################################################################
// Function Execute
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::JOB::
Execute(CACHED_ELIMINATION_MATRIX<T>* cem)
{
    auto& vl=cem->vector_list;
    int raw_op=op&raw_op_mask;
    MATRIX_MXN<T>* MO=0,*MA[3]={};
    //this initializes 
    for(int l=0;l<3;l++)
        if(arg_type[raw_op][l]==1 && a[l]>1)
            MA[l]=&cem->matrix_cache.Use(a[l]&raw_mask);
    if(arg_type[raw_op][3]==1 && o>1)
        MO=&cem->matrix_cache.Use(o);

    switch(raw_op)
    {
        case op_nop: break;
        case op_mat_inv:
            *MO=*MA[0];
            if(op&pseudo_inv) Pseudo_Inverse(*MO);
            else Inverse(*MO);
            break;
        case op_mat_mul:
            {
                int s0=0,s1=((a[1]^a[2])&use_neg)?-1:1;
                if(a[0]>=0) s0=(a[0]&use_neg)?-1:1;
                if(a[0]<0) MO->Resize(a[1]&use_trans?MA[1]->n:MA[1]->m,a[2]&use_trans?MA[2]->m:MA[2]->n);
                else if((a[0]&raw_mask)!=o) *MO=*MA[0];
                Times_MM(*MO,s0,*MA[1],a[1]&use_trans,*MA[2],a[2]&use_trans,s1);
            }
            break;
        case op_mat_add:
            {
                if(a[1]&use_trans)
                {
                    if(a[1]&use_neg) *MO=*MA[0]-MA[1]->Transposed();
                    else *MO=*MA[0]+MA[1]->Transposed();
                }
                else
                {
                    if(a[1]&use_neg) *MO=*MA[0]-*MA[1];
                    else *MO=*MA[0]+*MA[1];
                }
            }
            break;
        case op_vec_mul:
            {
                int s0=0,s1=((a[1]^a[2])&use_neg)?-1:1;
                if(a[0]>=0) s0=(a[0]&use_neg)?-1:1;
                auto& M=*MA[1];
                if(a[0]<0) vl(o).Resize(a[1]&use_trans?M.n:M.m);
                else if(a[0]!=o) vl(o)=vl(a[0]&raw_mask);
                Times_MV(vl(o),s1,M,a[1]&use_trans,vl(a[2]&raw_mask),s0);
            }
            break;
        case op_vec_add:
            if(a[1]&use_neg) vl(o)=vl(a[0])-vl(a[1]);
            else vl(o)=vl(a[0])+vl(a[1]);
            break;
        default: PHYSBAM_FATAL_ERROR();
    }

    for(int l=0;l<3;l++)
        if(arg_type[raw_op][l]>0 && a[l]>1)
        {
            bool rem=!--cem->data_refs[arg_type[raw_op][l]-1][a[l]&raw_mask];
            if(arg_type[raw_op][l]==1) cem->matrix_cache.Release(a[l]&raw_mask,rem);
            else if(rem) vl(a[l]&raw_mask).Clean_Memory();
        }
    if(arg_type[raw_op][3]==1 && o>1)
        cem->matrix_cache.Release(o,false);
}
//#####################################################################
// Function Compute_Job_Deps
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Compute_Job_Deps()
{
    provider[0].Resize(block_list.m,init_all,-1);
    provider[1].Resize(vector_list.m,init_all,-1);
    user[0].Resize(block_list.m);
    user[1].Resize(vector_list.m);
    for(int i=0;i<jobs.m;i++)
    {
        auto& j=jobs(i);
        if(arg_type[j.op&raw_op_mask][3]>0)
            provider[arg_type[j.op&raw_op_mask][3]-1](j.o)=i;
        for(int l=0;l<3;l++)
            if(arg_type[j.op&raw_op_mask][l]>0 && j.a[l]>=0)
                user[arg_type[j.op&raw_op_mask][l]-1](j.a[l]&raw_mask).Append(i);
    }
}
//#####################################################################
// Function Eliminate_Trans
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Eliminate_Trans()
{
    Foreach_Single_User_Pair([this](int p,int u)
        {
            auto& j=jobs(p);
            auto& k=jobs(u);
            if((j.op&raw_op_mask)==op_mat_mul && (k.op&raw_op_mask)==op_mat_add && (k.a[1]&use_trans) && j.o==(k.a[1]&raw_mask))
            {
                std::swap(j.a[1],j.a[2]);
                j.a[1]=Transposed(j.a[1]);
                j.a[2]=Transposed(j.a[2]);
                k.a[1]=Transposed(k.a[1]);
            }
        });
}
//#####################################################################
// Function Compute_Job_Deps
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Simplify_Jobs()
{
    for(auto& j:jobs)
    {
        switch(j.op&raw_op_mask)
        {
            case op_mat_mul:
            case op_vec_mul:
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
    Foreach_Single_User_Pair([this](int p,int u)
        {
            auto j=jobs(p);
            auto k=jobs(u);
            if((j.op&raw_op_mask)==op_mat_mul && (k.op&raw_op_mask)==op_mat_add && j.a[0]<0 && !(k.a[1]&use_trans))
            {
                if((k.a[1]&raw_mask)==j.o) std::swap(k.a[0],k.a[1]);
                PHYSBAM_ASSERT((k.a[0]&raw_mask)==j.o);
                j.a[0]=k.a[1];
                j.o=k.o;
                if(k.a[0]&use_neg) j.a[2]=Negate(j.a[2]);
                Remove_Job_Deps(p);
                Remove_Job_Deps(u);
                jobs(u)=j;
                jobs(p)={op_nop,{-1,-1,-1},-1};
                Add_Job_Deps(u);
            }
        });
}
//#####################################################################
// Function Provides_Argument
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Provides_Argument(const JOB& j,int a) const
{
    if(arg_type[j.op&raw_op_mask][a]==0) return -1;
    return provider[arg_type[j.op&raw_op_mask][a]-1](j.a[a]&raw_mask);
}
//#####################################################################
// Function Uses_Output
//#####################################################################
template<class T> const ARRAY<int>& CACHED_ELIMINATION_MATRIX<T>::
Uses_Output(const JOB& j) const
{
    static const ARRAY<int> empty;
    if(arg_type[j.op&raw_op_mask][3]==0) return empty;
    return user[arg_type[j.op&raw_op_mask][3]-1](j.o);
}
//#####################################################################
// Function Remove_Job_Deps
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Remove_Job_Deps(int j)
{
    const auto& k=jobs(j);
    if(arg_type[k.op&raw_op_mask][3]>0)
        provider[arg_type[k.op&raw_op_mask][3]-1](k.o)=-1;
    for(int i=0;i<3;i++)
        if(arg_type[k.op&raw_op_mask][i]>0 && k.a[i]>=0){
            auto& a=user[arg_type[k.op&raw_op_mask][i]-1](k.a[i]&raw_mask);
            int m=a.Find(j);
            PHYSBAM_ASSERT(m>=0);
            a.Remove_Index_Lazy(m);}
}
//#####################################################################
// Function Add_Job_Deps
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Add_Job_Deps(int j)
{
    const auto& k=jobs(j);
    if(arg_type[k.op&raw_op_mask][3]>0)
        provider[arg_type[k.op&raw_op_mask][3]-1](k.o)=j;
    for(int i=0;i<3;i++)
        if(arg_type[k.op&raw_op_mask][i]>0)
            user[arg_type[k.op&raw_op_mask][i]-1](k.a[i]&raw_mask).Append(j);
}
//#####################################################################
// Function Relabel
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Relabel()
{
    ARRAY<int> new_data[2],new_job(jobs.m,use_init,-1);
    int next[2]={num_orig_blocks,num_orig_vectors},next_job=0;
    for(int i=0;i<2;i++)
    {
        new_data[i].Resize(provider[i].m,use_init,-1);
        for(int j=0;j<next[i];j++)
            new_data[i](j)=j;
    }
    for(int i=0;i<jobs.m;i++)
    {
        auto& j=jobs(i);
        if((j.op&raw_op_mask)!=op_nop) new_job(i)=next_job++;
    }

    for(int t=0;t<2;t++)
    {
        for(int i=next[t];i<provider[t].m;i++)
        {
            int p=provider[t](i);
            if(p<0) continue;
            auto& j=jobs(p);
            if((j.op&raw_op_mask)!=op_nop && (j.op&raw_op_mask)!=op_mat_inv && j.a[0]>=0)
            {
                if(user[t](j.a[0]&raw_mask).m==1)
                {
                    new_data[t](i)=new_data[t](j.a[0]&raw_mask);
                    continue;
                }
            }
            new_data[t](i)=next[t]++;
        }
        data_refs[t]=new std::atomic<int>[next[t]];
        for(int i=0;i<next[t];i++) data_refs[t][i].store(0);
    }

    for(auto& a:rhs) if(a>=0) a=new_data[1](a);

    for(int i=0;i<jobs.m;i++)
    {
        if(new_job(i)<0) continue;
        JOB j=jobs(i);
        if(j.o>=0 && arg_type[j.op&raw_op_mask][3]>0)
            j.o=new_data[arg_type[j.op&raw_op_mask][3]-1](j.o);
        for(int i=0;i<3;i++)
            if(j.a[i]>=0 && arg_type[j.op&raw_op_mask][i]>0)
                j.a[i]=(j.a[i]&~raw_mask)|new_data[arg_type[j.op&raw_op_mask][i]-1](j.a[i]&raw_mask);
        jobs(new_job(i))=j;
    }
    jobs.Resize(next_job);

    for(int t=0;t<2;t++)
    {
        for(int i=0;i<provider[t].m;i++)
        {
            int p=provider[t](i);
            if(p<0) continue;
            p=new_job(p);
            for(auto u:user[t](i))
                dep_list.Append({p,new_job(u)});
        }
        provider[t].Clean_Memory();
        user[t].Clean_Memory();
    }

    for(int i=0;i<jobs.m;i++)
    {
        auto& j=jobs(i);
        for(int l=0;l<3;l++)
            if(j.a[l]>=0 && arg_type[j.op&raw_op_mask][l]>0)
                data_refs[arg_type[j.op&raw_op_mask][l]-1][j.a[l]&raw_mask]++;
    }
    for(auto& a:rhs) if(a>=0) data_refs[1][a]++;
}
//#####################################################################
// Function Register_Matrix_Block
//#####################################################################
template<class T> int CACHED_ELIMINATION_MATRIX<T>::
Create_Matrix_Block(bool sym)
{
    matrix_cache.entries.Add_End();
    return block_list.Append({sym,{block_list.m}});
}
//#####################################################################
// Function Add_Block_Matrix_Entry
//#####################################################################
template<class T> void CACHED_ELIMINATION_MATRIX<T>::
Add_Block_Matrix_Entry(int r,int c,int matrix_id)
{
    rows(r).Append({c,matrix_id});
    if(r!=c) rows(c).Append({r,Transposed(matrix_id)});
}
template struct CACHED_ELIMINATION_MATRIX<double>;
}