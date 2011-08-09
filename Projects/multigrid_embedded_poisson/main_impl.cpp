//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <boost/blank.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/variant/apply_visitor.hpp>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "Parsing/Parse_Command_Line.h"

#include "Build_And_Solve_Dirichlet_System.h"
#include "Build_And_Solve_Interface_System.h"
#include "Build_And_Solve_Neumann_System.h"
#include "Eval_Phi_Over_Fine_Grid.h"
#include "Lagrangify_Level_Set.h"
#include "RAND_MT19937_UNIFORM_REAL.h"

#ifndef PHYSBAM_NO_PETSC
#include <Jeffrey_Utilities/Petsc/CALL_AND_CHKERRQ.h>
#include <Jeffrey_Utilities/Petsc/SCOPED_FINALIZE.h>
#include <petsc.h>
#endif // #ifdef PHYSBAM_NO_PETSC

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace
{

template< class T, int D >
struct RUN_EXAMPLE_VISITOR;

} // namespace

template< class T, int D >
int main_impl(int argc, char* argv[])
{
    int main_result = 0;

    BASIC_TIMER timer;
    MAIN_PARAMS<T,D> main_params;

    std::cout << "Parsing command line..." << std::endl;
    timer.Restart();
    main_result = Parse_Command_Line(argc, argv, main_params);
    std::cout << "[Parsing command line...] " << timer.Elapsed() << " s" << std::endl;
    if(main_result != 0)
        return main_result;

#ifndef PHYSBAM_NO_PETSC
    const bool use_petsc =
        main_params.solver.solver_id == SOLVER_PARAMS::SOLVER_ID_PETSC_CG ||
        main_params.solver.solver_id == SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES;
    if(use_petsc) {
        std::cout << "Initializing PETSc...";
        std::cout.flush();
        timer.Restart();
        int petsc_argc = 0;
        char** petsc_argv = 0;
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( PetscInitialize(&petsc_argc, &petsc_argv, PETSC_NULL, PETSC_NULL) );
        std::cout << timer.Elapsed() << " s" << std::endl;
    }
    PHYSBAM_PETSC_SCOPED_FINALIZE_IF( use_petsc );
#endif // #ifndef PHYSBAM_NO_PETSC

    std::cout << "Initializing random number generator with seed "
              << main_params.general.rng_seed
              << "...";
    timer.Restart();
    boost::mt19937 rand_gen(main_params.general.rng_seed);
    static const boost::uniform_real<T> rand_dist(0, 1);
    typename RAND_MT19937_UNIFORM_REAL<T>::type rand(rand_gen, rand_dist);
    std::cout << timer.Elapsed() << " s" <<  std::endl;

    std::cout << "Evaluating level set function over (fine) grid...";
    std::cout.flush();
    timer.Restart();
    ARRAY<T> phi_of_fine_index((2 * As_Vector<int>(main_params.grid.n_cell) + 1).Product(), false); // uninit'ed
    Eval_Phi_Over_Fine_Grid(main_params, As_Array_View(phi_of_fine_index));
    std::cout << timer.Elapsed() << " s" << std::endl;

    main_result = boost::apply_visitor(
        RUN_EXAMPLE_VISITOR<T,D>(main_params, rand, phi_of_fine_index),
        main_params.example.problem
    );

    return main_result;
}

namespace
{

template< class T, int D >
struct RUN_EXAMPLE_VISITOR
{
    const MAIN_PARAMS<T,D>& main_params;
    typename RAND_MT19937_UNIFORM_REAL<T>::type& rand;
    const ARRAY_VIEW<const T> phi_of_fine_index;

    RUN_EXAMPLE_VISITOR(
        const MAIN_PARAMS<T,D>& main_params_,
        typename RAND_MT19937_UNIFORM_REAL<T>::type& rand_,
        const ARRAY_VIEW<const T> phi_of_fine_index_)
        : main_params(main_params_),
          rand(rand_),
          phi_of_fine_index(phi_of_fine_index_)
    { }

    typedef int result_type;

    int operator()(boost::blank) const
    { return Lagrangify_Level_Set(main_params, rand, phi_of_fine_index); }

    typedef typename EXAMPLE_PARAMS<T,D>::NEUMANN_PARAMS NEUMANN_PARAMS_TYPE;
    typedef typename EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS DIRICHLET_PARAMS_TYPE;
    typedef typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS INTERFACE_PARAMS_TYPE;

    int operator()(const NEUMANN_PARAMS_TYPE& problem) const
    { return Build_And_Solve_Neumann_System(problem, main_params, rand, phi_of_fine_index); }
    int operator()(const DIRICHLET_PARAMS_TYPE& problem) const
    { return Build_And_Solve_Dirichlet_System(problem, main_params, rand, phi_of_fine_index); }
    int operator()(const INTERFACE_PARAMS_TYPE& problem) const
    { return Build_And_Solve_Interface_System(problem, main_params, rand, phi_of_fine_index); }
};

} // namespace

#define EXPLICIT_INSTANTIATION( T, D ) \
template int main_impl<T,D>(int argc, char* argv[]);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
