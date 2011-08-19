//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_TWO_PHASE_FLOW_2D_CHORIN_PROJECT_HPP
#define PHYSBAM_PROJECTS_TWO_PHASE_FLOW_2D_CHORIN_PROJECT_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Functional/CONSTANT_FUNCTION.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/Multi_Index_X.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/RESULT_OF.h>
#include <Jeffrey_Utilities/SOLVER_PARAMS.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "../multigrid_embedded_poisson/Project_To_Zero_Divergence_With_Interface.h"

namespace PhysBAM
{

namespace Two_Phase_Flow_2D_Test
{

namespace Detail_Chorin_Project
{

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_KAPPA_OF_X_OF_CELL_INDEX,
    class T_MAC_VECTOR_FIELD
>
struct JUMP_P_OF_X_OF_CELL_INDEX;

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_KAPPA_OF_X,
    class T_MAC_VECTOR_FIELD
>
struct JUMP_P_OF_X;

} // namespace Detail_Chorin_Project

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_KAPPA_OF_X_OF_CELL_INDEX,
    class T_MAC_VECTOR_FIELD
>
void Chorin_Project(
    const unsigned int n_thread,
    const VECTOR<T,D>& min_x, const VECTOR<T,D>& max_x,
    const T dt,
    const T sigma, // coefficient of surface tension
    const T rho_negative, const T rho_positive, // density
    const T jump_mu, // viscosity
    const MULTI_INDEX_BOUND<D>& mac_cell_multi_index_bound,
    const T_PHI_OF_FINE_INDEX& phi_prev_of_fine_index,
    const T_PHI_OF_FINE_INDEX& phi_next_of_fine_index,
    const T_KAPPA_OF_X_OF_CELL_INDEX& kappa_of_x_of_cell_index, // curvature
    const T_MAC_VECTOR_FIELD& mac_vector_field)
{
    typedef Detail_Chorin_Project::JUMP_P_OF_X_OF_CELL_INDEX<
        T, D,
        T_PHI_OF_FINE_INDEX,
        T_KAPPA_OF_X_OF_CELL_INDEX,
        T_MAC_VECTOR_FIELD
    > JUMP_P_OF_X_OF_CELL_INDEX_;
    SOLVER_PARAMS solver_params;
    solver_params.relative_tolerance = 1e-12f;
    solver_params.print_diagnostics = true;
    Multigrid_Embedded_Poisson::Project_To_Zero_Divergence_With_Interface(
        n_thread,
        solver_params,
        min_x, max_x, mac_cell_multi_index_bound,
        phi_next_of_fine_index,
        Make_Constant_Function(1/rho_negative),
        Make_Constant_Function(1/rho_positive),
        JUMP_P_OF_X_OF_CELL_INDEX_(
            min_x, max_x, dt,
            sigma, jump_mu,
            mac_cell_multi_index_bound,
            phi_prev_of_fine_index,
            kappa_of_x_of_cell_index,
            mac_vector_field
        ),
        Make_Constant_Function(CONSTANT_FUNCTION<T>(static_cast<T>(0))),
        mac_vector_field,
        std::cout
    );
}

namespace Detail_Chorin_Project
{

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_KAPPA_OF_X_OF_CELL_INDEX,
    class T_MAC_VECTOR_FIELD
>
struct JUMP_P_OF_X_OF_CELL_INDEX
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        JUMP_P_OF_X_OF_CELL_INDEX,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, min_x ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, max_x ))
        (( typename T const, dt ))
        (( typename T const, sigma ))
        (( typename T const, jump_mu ))
        (( typename MULTI_INDEX_BOUND<D> const, mac_cell_multi_index_bound ))
        (( typename T_PHI_OF_FINE_INDEX const, phi_prev_of_fine_index ))
        (( typename T_KAPPA_OF_X_OF_CELL_INDEX const, kappa_of_x_of_cell_index ))
        (( typename T_MAC_VECTOR_FIELD const, mac_vector_field ))
    )
public:
    typedef typename RESULT_OF<
        const T_KAPPA_OF_X_OF_CELL_INDEX ( VECTOR<int,D> )
    >::type KAPPA_OF_X_TYPE;
    typedef JUMP_P_OF_X<
        T, D,
        T_PHI_OF_FINE_INDEX,
        KAPPA_OF_X_TYPE,
        T_MAC_VECTOR_FIELD
    > result_type;
    result_type operator()(const VECTOR<int,D> cell_multi_index) const
    {
        const MULTI_INDEX_BOUND<D> cell_multi_index_bound = mac_cell_multi_index_bound - 1;
        const VECTOR<T,D> dx = (max_x - min_x) / cell_multi_index_bound.max_multi_index;
        const VECTOR<T,D> x0 = Multi_Index_X(min_x, max_x, 2 * cell_multi_index_bound + 1, 2 * cell_multi_index);
        return result_type(
            dx/2,
            dt, sigma, jump_mu,
            phi_prev_of_fine_index,
            kappa_of_x_of_cell_index(cell_multi_index),
            mac_vector_field,
            cell_multi_index, x0
        );
    }
};

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_KAPPA_OF_X,
    class T_MAC_VECTOR_FIELD
>
struct JUMP_P_OF_X
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        JUMP_P_OF_X,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, dx_over_2 ))
        (( typename T const, dt ))
        (( typename T const, sigma ))
        (( typename T const, jump_mu ))
        (( typename T_PHI_OF_FINE_INDEX const, phi_prev_of_fine_index ))
        (( typename T_KAPPA_OF_X const, kappa_of_x ))
        (( typename T_MAC_VECTOR_FIELD const, mac_vector_field ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<int,D> )) const, cell_multi_index ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<T,D> )) const, x0 ))
    )
public:
    typedef T result_type;
    T operator()(VECTOR<T,D> x) const
    {
        typedef VECTOR<int,D> MULTI_INDEX_TYPE;
        // TODO: Add viscosity terms when jump_mu != 0.
        assert(jump_mu == 0);
        return dt * sigma * kappa_of_x(x);
#if 0
        x = (x - x0) / dx_over_2;
        MULTI_INDEX_TYPE fine_cell_multi_index = 2 * cell_multi_index;
        for(int d = 1; d <= D; ++d)
            if(x[d] < 0)
                --fine_cell_multi_index[d];
        x = PhysBAM::abs(x);
        T kappa_x = 0;
        BOOST_FOREACH( const MULTI_INDEX_TYPE fine_multi_offset, (STATIC_MULTI_INDEX_CUBE<D,0,1>()) ) {
            const MULTI_INDEX_TYPE fine_multi_index = fine_cell_multi_index + fine_multi_offset;
            const T weight = ((1 - fine_multi_offset) * (static_cast<T>(1) - x) + fine_multi_offset * x).Product();
            kappa_x += weight * kappa_of_fine_index(fine_multi_index);
        }
        return dt * sigma * kappa_x;
#endif // #if 0
    }
};

} // namespace Detail_Chorin_Project

} // namespace Two_Phase_Flow_2D_Test

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_TWO_PHASE_FLOW_2D_CHORIN_PROJECT_HPP
