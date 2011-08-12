//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <string>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception/exception.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include "../Params/MAIN_PARAMS.h"

#include "Parse_Command_Line.h"
#include "Parse_String_Example.h"
#include "Parse_String_Grid_BC.h"
#include "Parse_String_Grid_Box.h"
#include "Parse_String_Level_Set.h"
#include "Parse_String_N_Cell.h"
#include "Parse_String_Solver.h"
#include "Print_Help.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
int Parse_Command_Line(
    int argc, char* argv[],
    MAIN_PARAMS<T,D>& main_params)
{
    namespace bpo = boost::program_options;

    std::string grid_box_str;
    std::string n_cell_str;

    std::string level_set_str;

    std::string example_str;

    std::string grid_bc_str;

    std::string solver_str;

    bpo::options_description opts_desc("Options");
    opts_desc.add_options()
        ("help"           , "print this help message")
        ("dimension"      , bpo::value< unsigned int >()->
                            default_value(3)
                          , "spatial dimension" )
        ("scalar"         , bpo::value< std::string >()->
                            default_value("double")
                          , "scalar type (\"float\" | \"double\") (note that PETSc uses its own scalar type, specified when built)")
        ("seed"           , bpo::value< unsigned int >(&main_params.general.rng_seed)->
                            default_value(main_params.general.rng_seed)
                          , "random number seed to mt19937")
        ("n-thread"       , bpo::value< unsigned int >(&main_params.general.n_thread)->
                            default_value(main_params.general.n_thread)
                          , "number of threads to utilize")
        ("grid-box"       , bpo::value< std::string >(&grid_box_str)->
                            default_value("[-1,+1]^D")
                          , "grid domain box (\"[f,f]^D\" | \"[f,f] x [f,f] {x [f,f]}\")")
        ("n-cell"         , bpo::value< std::string >(&n_cell_str)
                          , "grid cell resolution (\"i^D\" | \"i x i {x i}\")")
        ("level-set"      ,  bpo::value< std::string >(&level_set_str)
                          , "level set specification (see below) for the boundary or interface surface")
        ("min-dist-to-vertex"
                          , bpo::value< float >(&main_params.level_set.min_dist_to_vertex)->
                            default_value(main_params.level_set.min_dist_to_vertex)
                          , "limits how close the level set rasterization gets to a grid vertex")
        ("example"        , bpo::value< std::string >(&example_str)
                          , "example specification (see below) defining the problem type, u, and beta")
        ("grid-bc"        , bpo::value< std::string >(&grid_bc_str)->
                            default_value("dirichlet")
                          , "grid boundary condition specification (\"dirichlet\" | \"neumann-offset\")")
        ("solver"         , bpo::value< std::string >(&solver_str)
                          , "solver specification (see below)")
        ("max-iterations" , bpo::value< unsigned int >(&main_params.solver.max_iterations)->
                            default_value(main_params.solver.max_iterations)
                          , "maximum number of solve iterations")
        ("rel-tol"        , bpo::value< float >(&main_params.solver.relative_tolerance)->
                            default_value(main_params.solver.relative_tolerance)
                          , "relative tolerance of solver")
        ("abs-tol"        , bpo::value< float >(&main_params.solver.absolute_tolerance)->
                            default_value(main_params.solver.absolute_tolerance)
                          , "absolute tolerance of solver")
        ("print-residuals", bpo::value< bool >(&main_params.solver.print_residuals)->
                            default_value(main_params.solver.print_residuals)
                          , "specifies whether to print the residual norm of the approximate solution at each solver iteration")
        ("precondition"   , bpo::value< bool >(&main_params.solver.precondition)->
                            default_value(main_params.solver.precondition)
                          , "specifies whether to use a preconditioner (only affects petsc-cg solves)")
        ("randomized-check"
                          , bpo::value< bool >(&main_params.general.randomized_check)->
                            default_value(main_params.general.randomized_check)
                          , "specifies whether to verify correctness of some operations via randomized input")
        ("embedded-surface-filename"
                          , bpo::value< std::string >(&main_params.output.embedded_surface_filename)
                          , "(output) filename for the Lagrangified embedded surface "
                            "(\"*.tri[.gz]\" | \"*.curve2d[.gz]\" | \"*.vtk\")")
        ("cell-local-embedded-surface-filename-format"
                          , bpo::value< std::string >(&main_params.output.cell_local_embedded_surface_filename_format)
                          , "(output) filename format (compatible with Boost.Format) for the cell-local Lagrangified embedded surface "
                            "(\"*.tri[.gz]\" | \"*.curve2d[.gz]\" | \"*.vtk\")")
        ("infty-norm-error-filename"
                          , bpo::value< std::string >(&main_params.output.infty_norm_error_filename)
                          , "(output) filename to write the infinity-norm errors");
    bpo::variables_map var_map;
    try {
        bpo::store(bpo::parse_command_line(argc, argv, opts_desc), var_map);
    }
    catch(boost::exception& e) {
        std::cerr << "ERROR: Exception thrown while parsing command line:\n"
                  << boost::diagnostic_information(e)
                  << "(Run with --help for usage.)"
                  << std::endl;
        return 1;
    }
    bpo::notify(var_map);
    if(var_map.count("help")) {
        std::cerr << opts_desc << std::endl;
        Print_Help(std::cerr);
        return 1;
    }

#define PARSE_STRING_AND_CHECK_RESULT( Name, Var, Expr ) \
    do { \
        const int result = Expr ; \
        if(result != 0) { \
            std::cerr << "ERROR: Unable to parse string parameter --" Name " : " << Var << '\n' \
                      << "(Run with --help for usage.)" \
                      << std::endl; \
            return result; \
        } \
    } while(false)

    PARSE_STRING_AND_CHECK_RESULT(
        "grid-box", grid_box_str,
        (Parse_String_Grid_Box(grid_box_str, main_params.grid.min_x, main_params.grid.max_x))
    );
    PARSE_STRING_AND_CHECK_RESULT(
        "n-cell", n_cell_str,
        (Parse_String_N_Cell(n_cell_str, main_params.grid.n_cell))
    );
    PARSE_STRING_AND_CHECK_RESULT(
        "level-set", level_set_str,
        (Parse_String_Level_Set(level_set_str, main_params.level_set.phi))
    );
    if(!example_str.empty()) {
        PARSE_STRING_AND_CHECK_RESULT(
            "example", example_str,
            (Parse_String_Example<T,D>(example_str, main_params.example.problem))
        );
        if(!main_params.example.Verify()) {
            std::cerr << "ERROR: Invalid u_id and/or beta_id in string parameter --example : " << example_str << '\n'
                      << "(Run with --help for usage.)"
                      << std::endl;
            return 1;
        }
    }
    PARSE_STRING_AND_CHECK_RESULT(
        "grid-bc", grid_bc_str,
        (Parse_String_Grid_BC(grid_bc_str, main_params.example.grid_bc_id))
    );
    if(!solver_str.empty())
        PARSE_STRING_AND_CHECK_RESULT(
            "solver", solver_str,
            (Parse_String_Solver(solver_str, main_params.solver))
        );

#undef PARSE_STRING_AND_CHECK_RESULT

    return 0;
}

#define EXPLICIT_INSTANTIATION( T, D ) \
template int Parse_Command_Line<T,D>( \
    int argc, char* argv[], \
    MAIN_PARAMS<T,D>& main_params);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
