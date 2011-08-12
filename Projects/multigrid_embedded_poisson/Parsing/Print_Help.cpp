//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PRINT_HELP_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PRINT_HELP_HPP

#include <ostream>

#include <boost/spirit/home/phoenix/core/argument.hpp>
#include <boost/spirit/home/phoenix/core/reference.hpp>
#include <boost/spirit/home/phoenix/operator/io.hpp>

#include <Jeffrey_Utilities/Parsing/BASIC_LEVEL_SET_GRAMMAR.h>
#include <Jeffrey_Utilities/Parsing/TRANSFORMED_LEVEL_SET_GRAMMAR.h>

#include "../Main_Examples.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

void Print_Help(std::ostream& out)
{
    namespace phx = boost::phoenix;
    namespace phx_arg = phx::arg_names;

    out << "Available level sets (dimension = 2):\n";
    BASIC_LEVEL_SET_GRAMMAR_BASE<2>::Visit_Grammar_Elements(
        phx::ref(out) << "  " << phx_arg::_1 << '\n'
    );
    out << "Available level sets (dimension = 3):\n";
    BASIC_LEVEL_SET_GRAMMAR_BASE<3>::Visit_Grammar_Elements(
        phx::ref(out) << "  " << phx_arg::_1 << '\n'
    );
    out << "Available level set transformations:\n";
    TRANSFORMED_LEVEL_SET_GRAMMAR_BASE::Visit_Grammar_Elements(
        phx::ref(out) << "  " << phx_arg::_1 << '\n'
    );

    out << "Available example specifications:\n"
        << "  neumann(<u>, <beta>)\n"
        << "  dirichlet(<u>, <beta>, <constraint>)\n"
        << "  interface((<u->, <beta->), (<u+>, <beta+>), <constraint>)\n";

    out << "Each of <u>, <beta>, <u->, <u+>, <beta->, <beta+> take the form:\n"
        << "  {<*_id>} |\n"
        << "  a*{<*_id>} |\n"
        << "  {<*_id>} + b |\n"
        << "  a*{<*_id>} + b\n"
        << "where a,b are scalars and the braces (\"{\", \"}\") are taken literally.\n";

    out << "Available u_id's (dimension = 2):\n";
    Visit_U_Examples_Str<2>(phx::ref(out) << "  " << phx_arg::_1 << " : " << phx_arg::_2 << '\n');
    out << "Available beta_id's (dimension = 2):\n";
    Visit_Beta_Examples_Str<2>(phx::ref(out) << "  " << phx_arg::_1 << " : " << phx_arg::_2 << '\n');

    out << "Available u_id's (dimension = 3):\n";
    Visit_U_Examples_Str<3>(phx::ref(out) << "  " << phx_arg::_1 << " : " << phx_arg::_2 << '\n');
    out << "Available beta_id's (dimension = 3):\n";
    Visit_Beta_Examples_Str<3>(phx::ref(out) << "  " << phx_arg::_1 << " : " << phx_arg::_2 << '\n');

    out << "Available constraints:\n"
        << "  single-cell\n"
        << "  double-cell{(<min-relative-indy-weight>)}\n"
        << "  aggregate{(<min-relative-indy-weight>)}\n";

    out << "Available solver specifications:\n"
        << "  [petsc|physbam]-[cg|minres]\n"
        << "  [mg|mgpcg]({<smoother>,}{<n-level>,}{<mu>,}(<nu1>,<nu11>,<nu12>,<nu13>){,(<nu2>,<nu21>,<nu22>,<nu23>)})\n";

    out << std::endl;
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PRINT_HELP_HPP
