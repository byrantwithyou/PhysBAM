//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_EVALUATE_CELL_LOCAL_INTEGRATIONS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_EVALUATE_CELL_LOCAL_INTEGRATIONS_HPP

#include <cassert>
#include <cmath>

#include <algorithm>

#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Geometry/Divide_Cell2.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_MULTI_POWERS_OVER_BOUNDARY_FOR_POISSON_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Geometry/INTEGRATE_MULTI_POWERS_OVER_VOLUME_FOR_POISSON_POLYTOPE_VISITOR.h>
#include <Jeffrey_Utilities/Math/Equal_Relative_Tolerance.h>
#include <Jeffrey_Utilities/Math/Eval_All_Grad_Multilinear_Basis_Dot_Grad_Multilinear_Basis_From_Monomial_Evals.h>
#include <Jeffrey_Utilities/Math/Eval_All_Multilinear_Basis_From_Monomial_Evals.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <Jeffrey_Utilities/VISITOR_SEQUENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T, int D,
    class T_PHI_OF_FINE_INDEX,
    class T_POLYTOPE_VISITOR
>
void Evaluate_Cell_Local_Integrations(
    const VECTOR<T,D> dx,
    const int domain_sign,
    const T_PHI_OF_FINE_INDEX& phi_of_fine_index,
    const float min_dist_to_vertex,
    const int sign_of_zero,
    const T_POLYTOPE_VISITOR& polytope_visitor,
    const VECTOR<int,D> cell_multi_index,
    VECTOR< T, STATIC_POW_C<3,D>::value >& monomials_integrated_over_volume,
    VECTOR< T, STATIC_POW_C<2,D>::value >& monomials_integrated_over_boundary,
    VECTOR< T, (1 << D) >& multilinear_basis_integrated_over_volume,
    VECTOR< T, (1 << D) >& multilinear_basis_integrated_over_boundary,
    VECTOR< VECTOR< T, (1 << D) >, (1 << D) >& grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume)
{
    static const T eps = std::numeric_limits<T>::epsilon();
    static_cast<void>(eps);

    const VECTOR<T,D> dx_over_2 = dx / 2;

    const T dv = dx.Product();
    static_cast<void>(dv);

    Divide_Cell2(
        dx_over_2,
        cell_multi_index,
        phi_of_fine_index,
        Make_Visitor_Sequence(
            Make_Integrate_Multi_Powers_Over_Volume_For_Poisson_Polytope_Visitor(
                domain_sign,
                Make_Compose_Function(
                    Make_Array_Wrapper_Function(monomials_integrated_over_volume),
                    STATIC_MULTI_INDEX_CUBE<D,0,2>()
                )
            ),
            Make_Integrate_Multi_Powers_Over_Boundary_For_Poisson_Polytope_Visitor(
                domain_sign,
                Make_Compose_Function(
                    Make_Array_Wrapper_Function(monomials_integrated_over_boundary),
                    STATIC_MULTI_INDEX_CUBE<D,0,1>()
                )
            ),
            polytope_visitor
        ),
        min_dist_to_vertex,
        sign_of_zero
    );

    const T volume_measure = monomials_integrated_over_volume[1];
    const T boundary_measure = monomials_integrated_over_boundary[1];
    static_cast<void>(volume_measure);
    static_cast<void>(boundary_measure);
    assert(0 < volume_measure && volume_measure <= dv);
    assert(boundary_measure >= 0);

    multilinear_basis_integrated_over_volume =
        Eval_All_Multilinear_Basis_From_Monomial_Evals(
            dx_over_2,
            Make_Compose_Function(
                Make_Array_Wrapper_Function(monomials_integrated_over_volume),
                STATIC_MULTI_INDEX_CUBE<D,0,2>()
            )
        );
    multilinear_basis_integrated_over_boundary =
        Eval_All_Multilinear_Basis_From_Monomial_Evals(
            dx_over_2,
            Make_Compose_Function(
                Make_Array_Wrapper_Function(monomials_integrated_over_boundary),
                STATIC_MULTI_INDEX_CUBE<D,0,1>()
            )
        );
    grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume =
        Eval_All_Grad_Multilinear_Basis_Dot_Grad_Multilinear_Basis_From_Monomial_Evals(
            dx_over_2,
            Make_Compose_Function(
                Make_Array_Wrapper_Function(monomials_integrated_over_volume),
                STATIC_MULTI_INDEX_CUBE<D,0,2>()
            )
        );

#ifndef NDEBUG
    {
        T alt_volume_measure = 0;
        T alt_boundary_measure = 0;
        for(int i = 1; i <= (1 << D); ++i) {
            assert(multilinear_basis_integrated_over_volume[i] > 0);
            assert(multilinear_basis_integrated_over_boundary[i] >= 0);
            assert(grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume[i][i] > 0);
            alt_volume_measure += multilinear_basis_integrated_over_volume[i];
            alt_boundary_measure += multilinear_basis_integrated_over_boundary[i];
            T sum = 0;
            T sum_abs = 0;
            for(int j = 1; j <= (1 << D); ++j) {
                assert(Equal_Relative_Tolerance<0>(
                    grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume[i][j],
                    grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume[j][i]
                ));
                sum += grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume[i][j];
                sum_abs += std::abs(grad_multilinear_basis_dot_grad_multilinear_basis_integrated_over_volume[i][j]);
            }
            assert(std::abs(sum) < (1<<D) * eps * std::max(sum_abs, (dv/(dx*dx)).Sum()));
        }
        assert(std::abs(alt_volume_measure - volume_measure) < (1<<D) * eps * dv);
        assert(std::abs(alt_boundary_measure - boundary_measure) < (1<<D) * eps * (dv/dx).Sum());
    }
#endif // #ifndef NDEBUG
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_EVALUATE_CELL_LOCAL_INTEGRATIONS_HPP
