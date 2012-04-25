//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_INTERFACE_BLOCK_HELPER
//#####################################################################

#ifndef __SYSTEM_INTERFACE_BLOCK_HELPER__
#define __SYSTEM_INTERFACE_BLOCK_HELPER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>

namespace PhysBAM{

template<class TV>
class SYSTEM_INTERFACE_BLOCK_HELPER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    
    VECTOR<MATRIX_MXN<T>,2> data;
    CELL_DOMAIN_INTERFACE<TV> *cdi;        
    CELL_MANAGER<TV> *cm;
    ARRAY<int> flat_diff;
    
    template<int d> 
    void Initialize(const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MANAGER<TV>& cm_input,
        CELL_DOMAIN_INTERFACE<TV> &cdi_input)
    {
        cm=&cm_input;
        cdi=&cdi_input;
        
        for(int i=0;i<s.diced.m;i++){
            const typename BASIS_STENCIL_UNIFORM<TV,d>::DICED& diced=s.diced(i);
            for(RANGE_ITERATOR<TV::m> it(cdi->coarse_range);it.Valid();it.Next())
                flat_diff.Append(cdi->Flatten_Diff(it.index+diced.index_offset));}
        
        flat_diff.Sort();
        flat_diff.Prune_Duplicates();
        
        for(int s=0;s<2;s++) data[s].Resize(cdi->interface_elements,flat_diff.m);
    }
    
    void Mark_Active_Cells(T tol)
    {
        for(int s=0;s<2;s++)
            for(int l=0;l<data[s].m;l++)
                for(int k=0;k<data[s].n;k++)
                    if(abs(data[s](l,k))>tol)
                        cm->Set_Active(cdi->Get_Flat_Base(l)+flat_diff(k),s);
                    else data[s](l,k)=0;
    }

    void Build_Matrix(VECTOR<SPARSE_MATRIX_FLAT_MXN<T>,2>& matrix)
    {
        if(!cdi->periodic_bc) PHYSBAM_FATAL_ERROR();

        for(int s=0;s<2;s++){
            SPARSE_MATRIX_FLAT_MXN<T>& M=matrix(s);
            MATRIX_MXN<T>& d=data[s];
            ARRAY<int>& comp_n=cm->compressed[s];
            int m=d.m;
            int n=cm->dofs[s];
            
            M.Reset(n);
            M.offsets.Resize(m+1);
            M.m=m;
            
            for(int row=0;row<m;row++){
                ARRAY<SPARSE_MATRIX_ENTRY<T> > entries;
                for(int j=0;j<d.n;j++){
                    int value=d(row,j);
                    if(value){
                        int column=comp_n(cdi->Get_Flat_Base(row)+flat_diff(j));
                        M.offsets(row+1)++;
                        entries.Append(SPARSE_MATRIX_ENTRY<T>(column,value));}}
                if(cdi->Is_Boundary_Element(row)) entries.Sort();
                M.A.Append_Elements(entries);}

            for(int i=0;i<M.offsets.m-1;i++) M.offsets(i+1)+=M.offsets(i);}
    }
};
}
#endif
