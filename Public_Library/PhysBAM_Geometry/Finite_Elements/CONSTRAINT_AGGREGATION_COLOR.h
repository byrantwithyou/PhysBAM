//#####################################################################
// Copyright 2012, Russell Howes.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTRAINT_AGGREGATION_COLOR
//#####################################################################
#ifndef __CONSTRAINT_AGGREGATION_COLOR__
#define __CONSTRAINT_AGGREGATION_COLOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
//#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_VECTOR_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR.h>



namespace PhysBAM{

    template<class TV> class GRID;
    
template<class TV>
class CONSTRAINT_AGGREGATION_COLOR:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV> VECTOR_T;

public:

    const GRID<TV>& grid;
    const CELL_DOMAIN_INTERFACE_COLOR<TV>* cdi;
    CELL_MANAGER_COLOR<TV>* cm;
    BASIS_INTEGRATION_UNIFORM_COLOR<TV,2>& biu;
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& full_constraint_matrix,condensed_constraint_matrix;
    
    ARRAY<T>& full_constraint_rhs;
    ARRAY<T> agg_constraint_rhs;
    ARRAY<T> weights;
    ARRAY<T> indices;
    ARRAY<int> neighbor_node_offsets;
    ARRAY<ARRAY<int> > corresponding_singlewide_constraints;
    ARRAY<VECTOR<T,2> > idof_color_pairs;
    ARRAY<int> idofs,is_idof_combined;
    ARRAY<ARRAY<int> >is_idof,is_non_idof;
//    ARRAY<VECTOR<int,2> > non_idofs;

    
    CONSTRAINT_AGGREGATION_COLOR(const GRID<TV>& grid_input,const CELL_DOMAIN_INTERFACE_COLOR<TV>* cdi_input, CELL_MANAGER_COLOR<TV>* cm_input,BASIS_INTEGRATION_UNIFORM_COLOR<TV,2>& biu_input,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& full_constraint_matrix_input,ARRAY<T>& full_constraint_rhs_input):grid(grid_input),cdi(cdi_input),cm(cm_input),biu(biu_input),full_constraint_matrix(full_constraint_matrix_input),condensed_constraint_matrix(),full_constraint_rhs(full_constraint_rhs_input){}
    
    void Set_Neighbor_Node_Offsets(){
        //Janky, can I make this nicer?
        TV_INT powers_of_three;
        int j=1;for(int i=TV::m-1;i>=0;i--){powers_of_three(i)=j;j*=3;}
        neighbor_node_offsets.Resize(j);
        
        TV_INT padded_cells;for(int i=0;i<TV::m;i++)padded_cells(i)=2*cdi->padding;padded_cells+=grid.numbers_of_cells;
        TV_INT a;
        
        a(TV::m-1)=1;for (int i=TV::m-2;i>=0;i--)a(i)=a(i+1)*padded_cells(i+1);
        
        for(int k=0;k<j;k++){
            int cell_temp=k;
            TV_INT cell_index;for(int i=TV::m-1;i>=0;i--){cell_index(i)=cell_temp % 3;cell_temp-=cell_index(i);cell_temp/=3;}
            cell_index-=1;
            neighbor_node_offsets(k)=cell_index.Dot(a);
           // std::cout<<k<<" "<<j<<" "<<neighbor_node_offsets(k)<<cell_index<<powers_of_three<<" "<<a<<"//";
        }
    }

//#####################################################################
    void Aggregate_Constraints(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix_uu,const ARRAY<int,TV_INT>& phi_color);
    void Build_Condensed_Constraint_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& agg_constraint_matrix);
    void Build_Nullspace_And_Specific_Constraint_Solution(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& a_matrix,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& z,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& z_transpose,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& z_t_a,SPARSE_MATRIX_FLAT_MXN<T>& ztaz,VECTOR_T& specific_solution);
    ARRAY<T>& RHS_Constraint_Agg(){return agg_constraint_rhs;}
private:
    void Get_Neighboring_Cells_From_Padded_Node(int node,int cube_radius,ARRAY<int>& neighbors);
    int Distance_Between_Node_And_Cell(int node,int cell);
    TV_INT Index_Of_Uncompressed_Node(int node);
//#####################################################################


};
}
#endif
