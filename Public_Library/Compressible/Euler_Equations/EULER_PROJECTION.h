//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_PROJECTION
//#####################################################################
#ifndef __EULER_PROJECTION__
#define __EULER_PROJECTION__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Compressible/Euler_Equations/EULER.h>
namespace PhysBAM{

template<class TV>
class EULER_PROJECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef TV_INT INDEX;
    typedef VECTOR<T,TV::m+2> TV_U;
    typedef typename ARRAY<T,TV_INT>::template REBIND<TV_U>::TYPE T_ARRAYS_U;

public:
    EULER_PROJECTION()
    {}

    virtual ~EULER_PROJECTION()
    {}
    //#####################################################################

    static void Compute_One_Over_rho_c_squared(const GRID<TV>& grid,const T_ARRAYS_U& U_ghost,const EOS<T>* eos,ARRAY<T,TV_INT>& one_over_rho_c_squared)
    {
        for(CELL_ITERATOR<TV>  iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
            INDEX cell_index=iterator.Cell_Index();
            const TV_U& U=U_ghost(cell_index);
            T rho=U(0); T c_squared=eos->c_squared(rho,EULER<TV>::e(U));
            one_over_rho_c_squared(cell_index)=(T)1/(rho*c_squared);}
    }
};
}
#endif
