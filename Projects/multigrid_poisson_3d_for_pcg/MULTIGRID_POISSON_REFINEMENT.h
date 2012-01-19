//#####################################################################
// Copyright 2009, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIGRID_POISSON_REFINEMENT
//#####################################################################
#ifndef __MULTIGRID_POISSON_REFINEMENT__
#define __MULTIGRID_POISSON_REFINEMENT__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include "BOX_ITERATOR.h"
#include "MULTIGRID_POISSON.h"

#ifndef MGPCG_UNOPTIMIZED
#include "../multigrid_poisson_3d_optimized_kernels/Coarsened_Discretization/Coarsened_Discretization_Helper.h"
#endif
namespace PhysBAM{

template<class T,int d> class MULTIGRID_POISSON;

template<class T,int d>
class MULTIGRID_POISSON_REFINEMENT
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<int,d> TV_INT;

public:

    MULTIGRID_POISSON<T,d>& fine;
    MULTIGRID_POISSON<T,d>& coarse;

    explicit MULTIGRID_POISSON_REFINEMENT(MULTIGRID_POISSON<T,d>& fine_input,MULTIGRID_POISSON<T,d>& coarse_input)
        :fine(fine_input),coarse(coarse_input)
    {}

    static MULTIGRID_POISSON<T,d>* Coarsened_Discretization(const MULTIGRID_POISSON<T,d>& fine_discretization)
    {
	for(int v=0;v<d;v++) if(fine_discretization.n(v)%2!=0) PHYSBAM_FATAL_ERROR("Number of cells in fine grid must be a multiple of two");
	MULTIGRID_POISSON<T,d>& coarse_discretization=*new MULTIGRID_POISSON<T,d>(fine_discretization.n/2,fine_discretization.h*2,fine_discretization.number_of_threads);

	// do boundary separately
	for(BOUNDARY_ITERATOR<d> coarse_iterator(coarse_discretization.padded_domain);coarse_iterator.Valid();coarse_iterator.Next()){
	    const T_INDEX& coarse_index=coarse_iterator.Index();
	    T_INDEX min_fine_index=(coarse_index-1)*2;
	    T_INDEX max_fine_index=min_fine_index+1;
	    min_fine_index=T_INDEX::Componentwise_Max(min_fine_index,fine_discretization.padded_domain.min_corner);
	    max_fine_index=T_INDEX::Componentwise_Min(max_fine_index,fine_discretization.padded_domain.max_corner);
	    coarse_discretization.cell_type(coarse_index)=MULTIGRID_POISSON<T,d>::NEUMANN_CELL_TYPE;
	    for(BOX_ITERATOR<d> fine_iterator(RANGE<T_INDEX>(min_fine_index,max_fine_index));fine_iterator.Valid();fine_iterator.Next()){
		const T_INDEX& fine_index=fine_iterator.Index();
		if(fine_discretization.cell_type(fine_index)==MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE)
		    coarse_discretization.cell_type(coarse_index)=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;
	    }
	}

#ifndef MGPCG_UNOPTIMIZED
	Coarsened_Discretization_Helper<T> helper(fine_discretization.n(1),fine_discretization.n(2),fine_discretization.n(3),&fine_discretization.cell_type(1,1,1),&coarse_discretization.cell_type(1,1,1));
	helper.Run_Parallel(fine_discretization.number_of_threads);
#else

	for(BOX_ITERATOR<d> coarse_iterator(coarse_discretization.unpadded_domain);coarse_iterator.Valid();coarse_iterator.Next()){
	    const T_INDEX& coarse_index=coarse_iterator.Index();
	    const T_INDEX base_fine_index=(coarse_index-1)*2;
	    coarse_discretization.cell_type(coarse_index)=MULTIGRID_POISSON<T,d>::NEUMANN_CELL_TYPE;
	    for(BOX_ITERATOR<d> fine_iterator(RANGE<T_INDEX>(base_fine_index,base_fine_index+1));fine_iterator.Valid();fine_iterator.Next()){
		const T_INDEX& fine_index=fine_iterator.Index();
		if(fine_discretization.cell_type(fine_index)==MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE)
		    coarse_discretization.cell_type(coarse_index)=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;
		else if(fine_discretization.cell_type(fine_index)==MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE && 
		    coarse_discretization.cell_type(coarse_index)!=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE)
		    coarse_discretization.cell_type(coarse_index)=MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE;}}

#endif
	coarse_discretization.Initialize();
	return &coarse_discretization;
    }
//#####################################################################
    void Transfer_Residual_To_Coarse_Grid();
    void Transfer_Correction_To_Fine_Grid();
//#####################################################################
};
}
#endif
