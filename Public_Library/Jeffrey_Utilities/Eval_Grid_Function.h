//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_EVAL_GRID_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_EVAL_GRID_FUNCTION_HPP

#include <cassert>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Algorithm/For_Each.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_X.h>

namespace PhysBAM
{

//#####################################################################
//#####################################################################

template< class T, int D, class T_F_OF_X >
inline void
Eval_Grid_Function(
    const VECTOR<T,D>& min_x, const VECTOR<T,D>& max_x,
    const MULTI_INDEX_BOUND<D> multi_index_bound,
    T_F_OF_X f_of_x,
    ARRAY_VIEW<T> f_of_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    const int n = multi_index_bound.Size();
    assert(f_of_index.Size() == n);
    for(int linear_index = 1; linear_index <= n; ++linear_index) {
        const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
        const VECTOR<T,D> x = Multi_Index_X(min_x, max_x, multi_index_bound, multi_index);
        f_of_index(linear_index) = f_of_x(x);
    }
}

//#####################################################################
//#####################################################################

namespace Detail_Eval_Grid_Function
{

template< class T, int D, class T_F_OF_X >
struct EVAL_GRID_FUNCTION_HELPER;

} // namespace Detail_Eval_Grid_Function

template< class T, int D, class T_F_OF_X >
inline void
Eval_Grid_Function_MT(
    const unsigned int n_thread,
    const VECTOR<T,D>& min_x, const VECTOR<T,D>& max_x,
    const MULTI_INDEX_BOUND<D>& multi_index_bound,
    const T_F_OF_X& f_of_x,
    ARRAY_VIEW<T> f_of_index)
{
    typedef Detail_Eval_Grid_Function::EVAL_GRID_FUNCTION_HELPER<
        T, D, T_F_OF_X
    > EVAL_GRID_FUNCTION_HELPER_;
    assert(n_thread >= 1);
    assert(f_of_index.Size() == multi_index_bound.Size());
    For_Each_MT(
        n_thread,
        1, multi_index_bound.Size(),
        EVAL_GRID_FUNCTION_HELPER_(
            min_x, max_x, multi_index_bound,
            f_of_x, f_of_index
        )
    );
}

//#####################################################################
//#####################################################################

namespace Detail_Eval_Grid_Function
{

template< class T, int D, class T_F_OF_X >
struct EVAL_GRID_FUNCTION_HELPER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        EVAL_GRID_FUNCTION_HELPER,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, min_x ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, max_x ))
        (( typename MULTI_INDEX_BOUND<D> const, multi_index_bound ))
        (( typename T_F_OF_X const, f_of_x ))
        (( typename ARRAY_VIEW<T>&, f_of_index ))
    )
public:
    typedef void result_type;
    void operator()(const int linear_index) const
    {
        const VECTOR<int,D> multi_index = multi_index_bound.Multi_Index(linear_index);
        const VECTOR<T,D> x = Multi_Index_X(min_x, max_x, multi_index_bound, multi_index);
        f_of_index(linear_index) = f_of_x(x);
    }
};

} // namespace Detail_Eval_Grid_Function

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_EVAL_GRID_FUNCTION_HPP
