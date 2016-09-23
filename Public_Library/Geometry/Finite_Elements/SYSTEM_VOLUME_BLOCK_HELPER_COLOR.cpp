//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>
#include <map>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> template<int d0,class ...Args> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Initialize(CELL_DOMAIN_INTERFACE_COLOR<TV> &cdi_input,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,CELL_MANAGER_COLOR<TV>& cm0,Args&& ...args)
{
    cdi=&cdi_input;
    Store_Cell_Manager(s0,cm0,args...);

    ARRAY<int> diffs;
    for(int i=0;i<s0.diced.m;i++)
        Initialize_Helper(s0.diced(i).subcell,diffs,s0.diced(i).index_offset,args...);

    Prune_Duplicates(flat_diff);
    flat_diff.Sort(LEXICOGRAPHIC_COMPARE());
    
    data.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++) data(c).Resize(cdi->flat_size,flat_diff.m);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> template<int d1,class ...Args> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Initialize_Helper(int subcell,ARRAY<int>& diffs,const TV_INT& index_offset,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,CELL_MANAGER_COLOR<TV>& cm1,Args&& ...args)
{
    for(int j=0;j<s1.diced.m;j++)
        if(int new_subcell=subcell&s1.diced(j).subcell){
            diffs.Append(cdi->Flatten_Diff(s1.diced(j).index_offset-index_offset));
            Initialize_Helper(new_subcell,diffs,index_offset,args...);
            diffs.Pop();}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Initialize_Helper(int subcell,ARRAY<int>& diffs,const TV_INT& index_offset)
{
    flat_diff.Append(diffs);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> template<int d1,class ...Args> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Store_Cell_Manager(const BASIS_STENCIL_UNIFORM<TV,d1>& s1,CELL_MANAGER_COLOR<TV>& cm1,Args&& ...args)
{
    cm.Append(&cm1);
    Store_Cell_Manager(args...);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Store_Cell_Manager()
{
}
//#####################################################################
// Function Mark_Active_Cells
//#####################################################################
template<class TV> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Mark_Active_Cells(T tol)
{
    // add ghost rows to material ones 
    for(CELL_ITERATOR<TV> it(cdi->grid,cdi->padding,GRID<TV>::GHOST_REGION,-1);it.Valid();it.Next()){
        int i=cdi->Flatten(it.index);
        int r=cdi->remap(i);
        for(int c=0;c<cdi->colors;c++)
            for(int k=0;k<data(c).n;k++){
                data(c)(r,k)+=data(c)(i,k);
                data(c)(i,k)=0;}}
    for(int c=0;c<cdi->colors;c++)
        for(int l=0;l<data(c).m;l++)
            for(int k=0;k<data(c).n;k++)
                if(abs(data(c)(l,k))>tol){
                    cm(0)->Set_Active(l,c);
                    for(int p=0;p<flat_diff(k).m;p++)
                        cm(p+1)->Set_Active(l+flat_diff(k)(p),c);}
                else data(c)(l,k)=0;
}
//#####################################################################
// Function Build_Matrix
//#####################################################################
template<class TV> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Build_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix)
{
    matrix.Resize(cdi->colors);
    PHYSBAM_ASSERT(cm.m==2);

    for(int c=0;c<cdi->colors;c++){
        SPARSE_MATRIX_FLAT_MXN<T>& M=matrix(c);
        MATRIX_MXN<T>& d=data(c);
        ARRAY<int>& comp_m=cm(0)->compressed(c);
        ARRAY<int>& comp_n=cm(1)->compressed(c);
        int m=cm(0)->dofs(c);
        int n=cm(1)->dofs(c);
            
        M.Reset(n);
        M.offsets.Resize(m+1);
        M.m=m;
            
        for(int i=0;i<d.m;i++){
            if(cdi->Is_Outside_Cell(i)) continue;
            int row=comp_m(i);
            if(row>=0){
                ARRAY<SPARSE_MATRIX_ENTRY<T> > entries;
                for(int j=0;j<d.n;j++){
                    T value=d(i,j);
                    if(value){
                        int column=comp_n(i+flat_diff(j)(0));
                        M.offsets(row+1)++;
                        entries.Append(SPARSE_MATRIX_ENTRY<T>(column,value));}}
                if(cdi->Is_Boundary_Cell(i)) entries.Sort();
                M.A.Append_Elements(entries);}}

        for(int i=0;i<M.offsets.m-1;i++) M.offsets(i+1)+=M.offsets(i);}
}
//#####################################################################
// Function Build_Matrix
//#####################################################################
template<class TV> template<class ...Args> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Build_Matrix_With_Contract(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix,Args&& ...args)
{
    matrix.Resize(cdi->colors);
    ARRAY<ARRAY<T> > contract_row;
    ARRAY<ARRAY<VECTOR<int,2> > > instructions;
    ARRAY<ARRAY<const ARRAY<T>*> > vecs(cdi->colors);
    ARRAY<CELL_MANAGER_COLOR<TV>*> pruned_cm(cm);
    ARRAY<int> contract_index;
    ARRAY<ARRAY<int> > new_flat_diff;
    Setup_Contract_Instructions(contract_row,instructions,vecs,flat_diff,new_flat_diff,pruned_cm,contract_index,args...);

    ARRAY<ARRAY_VIEW<T> > contract_row_view(contract_row.m+1);
    for(int i=0;i<contract_row.m;i++)
        contract_row_view(i+1).Set(contract_row(i));

    for(int c=0;c<cdi->colors;c++){
        SPARSE_MATRIX_FLAT_MXN<T>& M=matrix(c);
        MATRIX_MXN<T>& d=data(c);
        ARRAY<int>* comp_m=&pruned_cm(0)->compressed(c);
        ARRAY<int>* comp_n=0;
        int m=pruned_cm(0)->dofs(c),n=-1;
        for(int i=1;i<pruned_cm.m;i++)
            if(pruned_cm(i)){
                comp_n=&pruned_cm(i)->compressed(c);
                n=pruned_cm(i)->dofs(c);}
    
        M.Reset(n);
        M.offsets.Resize(m+1);
        M.m=m;
            
        for(int i=0;i<d.m;i++){
            if(cdi->Is_Outside_Cell(i)) continue;
            int row=(*comp_m)(i);
            if(row>=0){
                contract_row_view(0).Set(d.x.Array_View(d.n*i,d.n));
                for(int j=0;j<instructions.m;j++){
                    ARRAY<int>& comp=cm(contract_index(j))->compressed(c);
                    contract_row_view(j+1).Fill(0);
                    for(int k=0;k<instructions(j).m;k++){
                        if(T value=contract_row_view(j)(k)){
                            int ind=comp(i+instructions(j)(k).x);
                            contract_row_view(j+1)(instructions(j)(k).y)+=
                                (*vecs(c)(j))(ind)*value;}}}
                ARRAY_VIEW<T> matrix_row(contract_row_view.Last());

                ARRAY<SPARSE_MATRIX_ENTRY<T> > entries;
                for(int j=0;j<matrix_row.m;j++){
                    T value=matrix_row(j);
                    if(value){
                        int column=(*comp_n)(i+new_flat_diff(j)(0));
                        M.offsets(row+1)++;
                        assert(column>=0);
                        entries.Append(SPARSE_MATRIX_ENTRY<T>(column,value));}}
                if(cdi->Is_Boundary_Cell(i)) entries.Sort();
                M.A.Append_Elements(entries);}}

        for(int i=0;i<M.offsets.m-1;i++) M.offsets(i+1)+=M.offsets(i);}
}
//#####################################################################
// Function Build_Matrix
//#####################################################################
template<class TV> template<class ...Args> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Setup_Contract_Instructions(ARRAY<ARRAY<T> >& contract_row,ARRAY<ARRAY<VECTOR<int,2> > >& instructions,
    ARRAY<ARRAY<const ARRAY<T>*> >& vecs,const ARRAY<ARRAY<int> >& flat_diff,ARRAY<ARRAY<int> >& new_flat_diff,
    ARRAY<CELL_MANAGER_COLOR<TV>*>& pruned_cm,ARRAY<int>& contract_index,int index,const ARRAY<ARRAY<T> >& vec,
    Args&& ...args)
{
    contract_index.Append(index);
    pruned_cm(index)=0;
    instructions.Append(ARRAY<VECTOR<int,2> >());
    instructions.Last().Resize(flat_diff.m);
    std::map<ARRAY<int>,int,LEXICOGRAPHIC_COMPARE> flat_diff_map;
    for(int i=0;i<flat_diff.m;i++){
        ARRAY<int> a=flat_diff(i);
        instructions.Last()(i).x=a(index-1);
        a(index-1)=0;
        instructions.Last()(i).y=flat_diff_map.emplace(a,flat_diff_map.size()).first->second;}
    contract_row.Append(ARRAY<T>());
    contract_row.Last().Resize(flat_diff_map.size());
    for(int c=0;c<cdi->colors;c++)
        vecs(c).Append(&vec(c));

    ARRAY<ARRAY<int> > flatter_diff;
    for(const auto& it:flat_diff_map) flatter_diff.Append(it.first);

    Setup_Contract_Instructions(contract_row,instructions,vecs,flatter_diff,new_flat_diff,pruned_cm,contract_index,args...);
}
//#####################################################################
// Function Build_Matrix
//#####################################################################
template<class TV> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Setup_Contract_Instructions(ARRAY<ARRAY<T> >& contract_row,ARRAY<ARRAY<VECTOR<int,2> > >& instructions,
    ARRAY<ARRAY<const ARRAY<T>*> >& vecs,const ARRAY<ARRAY<int> >& flat_diff,ARRAY<ARRAY<int> >& new_flat_diff,
    ARRAY<CELL_MANAGER_COLOR<TV>*>& pruned_cm,ARRAY<int>& contract_index)
{
    new_flat_diff=flat_diff;
}
namespace PhysBAM{
template class SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >;
template class SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >;
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >::Build_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<double>,int>&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >::Mark_Active_Cells(double);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >::Build_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<double>,int>&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >::Mark_Active_Cells(double);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,CELL_MANAGER_COLOR<VECTOR<double,2> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,CELL_MANAGER_COLOR<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,CELL_MANAGER_COLOR<VECTOR<double,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,CELL_MANAGER_COLOR<VECTOR<double,2> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,CELL_MANAGER_COLOR<VECTOR<double,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,CELL_MANAGER_COLOR<VECTOR<double,3> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,CELL_MANAGER_COLOR<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,CELL_MANAGER_COLOR<VECTOR<double,3> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,CELL_MANAGER_COLOR<VECTOR<double,3> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,CELL_MANAGER_COLOR<VECTOR<double,3> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,CELL_MANAGER_COLOR<VECTOR<float,2> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,CELL_MANAGER_COLOR<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,CELL_MANAGER_COLOR<VECTOR<float,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,CELL_MANAGER_COLOR<VECTOR<float,2> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,CELL_MANAGER_COLOR<VECTOR<float,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,CELL_MANAGER_COLOR<VECTOR<float,3> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,CELL_MANAGER_COLOR<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,CELL_MANAGER_COLOR<VECTOR<float,3> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,CELL_MANAGER_COLOR<VECTOR<float,3> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,CELL_MANAGER_COLOR<VECTOR<float,3> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >::Build_Matrix_With_Contract<int,ARRAY<ARRAY<double,int>,int>&>(ARRAY<SPARSE_MATRIX_FLAT_MXN<double>,int>&,int&&,ARRAY<ARRAY<double,int>,int>&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,CELL_MANAGER_COLOR<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,CELL_MANAGER_COLOR<VECTOR<double,2> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,CELL_MANAGER_COLOR<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,CELL_MANAGER_COLOR<VECTOR<double,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >::Build_Matrix_With_Contract<int,ARRAY<ARRAY<double,int>,int>&>(ARRAY<SPARSE_MATRIX_FLAT_MXN<double>,int>&,int&&,ARRAY<ARRAY<double,int>,int>&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,CELL_MANAGER_COLOR<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,CELL_MANAGER_COLOR<VECTOR<double,3> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,CELL_MANAGER_COLOR<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,CELL_MANAGER_COLOR<VECTOR<double,3> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >::Build_Matrix_With_Contract<int,ARRAY<ARRAY<float,int>,int>&>(ARRAY<SPARSE_MATRIX_FLAT_MXN<float>,int>&,int&&,ARRAY<ARRAY<float,int>,int>&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,CELL_MANAGER_COLOR<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,CELL_MANAGER_COLOR<VECTOR<float,2> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,CELL_MANAGER_COLOR<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,CELL_MANAGER_COLOR<VECTOR<float,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >::Build_Matrix_With_Contract<int,ARRAY<ARRAY<float,int>,int>&>(ARRAY<SPARSE_MATRIX_FLAT_MXN<float>,int>&,int&&,ARRAY<ARRAY<float,int>,int>&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,CELL_MANAGER_COLOR<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,CELL_MANAGER_COLOR<VECTOR<float,3> >&>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,CELL_MANAGER_COLOR<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,CELL_MANAGER_COLOR<VECTOR<float,3> >&);
}
