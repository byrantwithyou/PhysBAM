//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <string>

#include <boost/spirit/home/phoenix/core/reference.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/spirit/home/qi/action/action.hpp>
#include <boost/spirit/home/qi/char/char_class.hpp>
#include <boost/spirit/home/qi/operator/alternative.hpp>
#include <boost/spirit/home/qi/parse.hpp>
#include <boost/spirit/home/qi/string/lit.hpp>

#include "../Params/EXAMPLE_PARAMS.h"

#include "Parse_String_Grid_BC.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

int Parse_String_Grid_BC(
    const std::string& grid_bc_str,
    EXAMPLE_PARAMS_BASE::GRID_BC_ID& grid_bc_id)
{
    namespace phx = boost::phoenix;
    namespace qi  = boost::spirit::qi;

    std::string::const_iterator it = grid_bc_str.begin();
    const bool b = qi::phrase_parse(
        it, grid_bc_str.end(),
        qi::lit("neumann")  [phx::ref(grid_bc_id) = EXAMPLE_PARAMS_BASE::GRID_BC_ID_NEUMANN  ] |
        qi::lit("dirichlet")[phx::ref(grid_bc_id) = EXAMPLE_PARAMS_BASE::GRID_BC_ID_DIRICHLET],
        qi::standard::blank
    );
    if(b && it == grid_bc_str.end())
        return 0;

    return 1;
}

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM
