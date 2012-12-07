//#####################################################################
// Copyright 2012, Russell Howes.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTRAINT_AGGREGATION_COLOR
//#####################################################################
#ifndef __CONSTRAINT_AGGREGATION_COLOR__
#define __CONSTRAINT_AGGREGATION_COLOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM_COLOR.h>



namespace PhysBAM{

    template<class TV> class GRID;
    
template<class TV>
class CONSTRAINT_AGGREGATION_COLOR:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:

    const GRID<TV>& grid;
    const CELL_DOMAIN_INTERFACE_COLOR<TV>* cdi;
    CELL_MANAGER_COLOR<TV>* cm;
    BASIS_INTEGRATION_UNIFORM_COLOR<TV,2>& biu;
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& full_constraint_matrix;
    ARRAY<T>& full_constraint_rhs;
    ARRAY<T> weights;
    ARRAY<T> indices;
    ARRAY<ARRAY<int> > corresponding_singlewide_constraints;
    ARRAY<VECTOR<T,2> > idof_color_pairs;
    ARRAY<int> idofs;

    
    CONSTRAINT_AGGREGATION_COLOR(const GRID<TV>& grid_input,const CELL_DOMAIN_INTERFACE_COLOR<TV>* cdi_input, CELL_MANAGER_COLOR<TV>* cm_input,BASIS_INTEGRATION_UNIFORM_COLOR<TV,2>& biu_input,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& full_constraint_matrix_input,ARRAY<T>& full_constraint_rhs_input):grid(grid_input),cdi(cdi_input),cm(cm_input),biu(biu_input),full_constraint_matrix(full_constraint_matrix_input),full_constraint_rhs(full_constraint_rhs_input){}

//#####################################################################
    void Aggregate_Constraints(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix_uu,SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV>& helper_qu,const ARRAY<int,TV_INT>& phi_color);
    void Build_Condensed_Constraint_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& agg_constraint_matrix,ARRAY<T>& constraint_rhs);
private:
    void Get_Neighboring_Cells_From_Padded_Node(int node,int cube_radius,ARRAY<int>& neighbors);
    int Distance_Between_Node_And_Cell(int node,int cell);
    TV_INT Index_Of_Uncompressed_Node(int node);
//#####################################################################


};
}
#endif
