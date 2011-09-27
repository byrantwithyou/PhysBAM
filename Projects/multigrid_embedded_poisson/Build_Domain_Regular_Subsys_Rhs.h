//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_DOMAIN_REGULAR_SUBSYS_RHS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_DOMAIN_REGULAR_SUBSYS_RHS_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Grid/Visit_Cells_With_Sign_Via_Cell_Sign.h>
#include <Jeffrey_Utilities/Grid/VISIT_IF_SIGN_PREDICATE_GRID_VISITOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace Detail_Build_Domain_Regular_Subsys_Rhs
{

template< class T, int D, class T_F_OF_INDEX >
struct BUILD_DOMAIN_REGULAR_SUBSYS_RHS_VISITOR;

} // namespace Detail_Build_Domain_Regular_Subsys_Rhs

template< class T, int D, class T_SIGN_OF_CELL_INDEX, class T_F_OF_INDEX >
void Build_Domain_Regular_Subsys_Rhs(
    const unsigned int n_thread,
    const T dv,
    const MULTI_INDEX_BOUND<D>& cell_multi_index_bound,
    const int domain_sign,
    const T_SIGN_OF_CELL_INDEX& sign_of_cell_index,
    const T_F_OF_INDEX& f_of_index,
    ARRAY_VIEW<T> system_rhs)
{
    // TODO: Need to stripe the domain...
    typedef Detail_Build_Domain_Regular_Subsys_Rhs::BUILD_DOMAIN_REGULAR_SUBSYS_RHS_VISITOR<
        T, D, T_F_OF_INDEX
    > BUILD_DOMAIN_REGULAR_SUBSYS_RHS_VISITOR_;
    assert(n_thread >= 1);
    assert(system_rhs.Size() >= (cell_multi_index_bound + 1).Size());
    Visit_Cells_With_Sign_Via_Cell_Sign(
        cell_multi_index_bound.Size(),
        sign_of_cell_index,
        Make_Visit_If_Sign_Predicate_Grid_Visitor(
            Make_Equal_Function(domain_sign),
            BUILD_DOMAIN_REGULAR_SUBSYS_RHS_VISITOR_(
                dv, cell_multi_index_bound,
                f_of_index,
                system_rhs
            )
        )
    );
}

namespace Detail_Build_Domain_Regular_Subsys_Rhs
{

template< class T, int D, class T_F_OF_INDEX >
struct BUILD_DOMAIN_REGULAR_SUBSYS_RHS_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        BUILD_DOMAIN_REGULAR_SUBSYS_RHS_VISITOR,
        (( typename T const, dv ))
        (( typename MULTI_INDEX_BOUND<D> const, cell_multi_index_bound ))
        (( typename T_F_OF_INDEX const, f_of_index ))
        (( typename ARRAY_VIEW<T>&, system_rhs ))
    )
public:
    typedef void result_type;
    void operator()(const int cell_linear_index) const
    {
        typedef VECTOR<int,D> MULTI_INDEX_TYPE;
        const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
        const MULTI_INDEX_TYPE cell_multi_index = cell_multi_index_bound.Multi_Index(cell_linear_index);
        BOOST_FOREACH(
            const MULTI_INDEX_TYPE multi_index,
            (MULTI_INDEX_CUBE<D,0,1>(cell_multi_index))
        ) {
            const int linear_index = multi_index_bound.Linear_Index(multi_index);
            system_rhs(linear_index) += f_of_index(multi_index) * (dv / (1 << D));
        }
    }
};

} // namespace Detail_Build_Domain_Regular_Subsys_Rhs

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_DOMAIN_REGULAR_SUBSYS_RHS_HPP
