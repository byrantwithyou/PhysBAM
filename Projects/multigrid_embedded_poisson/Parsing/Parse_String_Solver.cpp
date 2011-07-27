//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <algorithm>
#include <string>

#include <boost/fusion/container/vector/vector10.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/spirit/home/phoenix/core/reference.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/spirit/home/phoenix/statement/sequence.hpp>
#include <boost/spirit/home/qi/action/action.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/char/char_class.hpp>
#include <boost/spirit/home/qi/nonterminal/rule.hpp>
#include <boost/spirit/home/qi/numeric/uint.hpp>
#include <boost/spirit/home/qi/operator/alternative.hpp>
#include <boost/spirit/home/qi/operator/optional.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/parse.hpp>
#include <boost/spirit/home/qi/string/lit.hpp>
#include <boost/spirit/home/support/argument.hpp>
#include <boost/spirit/home/support/context.hpp>

#include "../Params/SOLVER_PARAMS.h"

#include "Parse_String_Solver.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

int Parse_String_Solver(
    const std::string& solver_str,
    SOLVER_PARAMS& solver)
{
    namespace phx = boost::phoenix;
    namespace qi  = boost::spirit::qi;

    static const qi::rule<
        std::string::const_iterator,
        SOLVER_PARAMS::MULTIGRID_PARAMS::SMOOTHER_ID ( ),
        qi::standard::blank_type
    > smoother =
        (qi::lit("gs") | qi::lit("gauss-seidel"))
        [ qi::_val = SOLVER_PARAMS::MULTIGRID_PARAMS::SMOOTHER_ID_GAUSS_SEIDEL ]
      | (qi::lit("jacobi") | qi::lit("weighted-jacobi"))
        [ qi::_val = SOLVER_PARAMS::MULTIGRID_PARAMS::SMOOTHER_ID_WEIGHTED_JACOBI ];

    typedef boost::fusion::vector4<
        unsigned int, unsigned int, unsigned int, unsigned int
    > uint4_type;
    static const qi::rule<
        std::string::const_iterator,
        uint4_type ( ),
        qi::standard::blank_type
    > uint4 =
        '('
            >> qi::uint_ >> ','
            >> qi::uint_ >> ','
            >> qi::uint_ >> ','
            >> qi::uint_ >>
        ')';

    SOLVER_PARAMS::SOLVER_ID& solver_id = solver.solver_id;
    unsigned int& n_level = (solver.multigrid.n_level = 0);
    unsigned int& mu = (solver.multigrid.mu = 1);
    SOLVER_PARAMS::MULTIGRID_PARAMS::SMOOTHER_ID& smoother_id =
        (solver.multigrid.smoother_id = SOLVER_PARAMS::MULTIGRID_PARAMS::SMOOTHER_ID_GAUSS_SEIDEL);
    typedef boost::fusion::vector4<
        unsigned int&, unsigned int&, unsigned int&, unsigned int&
    > uintref4_type;
    uintref4_type nu_pre_restrict(
        solver.multigrid.nu_pre_restrict.n,
        solver.multigrid.nu_pre_restrict.n_pre_full_embedded,
        solver.multigrid.nu_pre_restrict.n_full,
        solver.multigrid.nu_pre_restrict.n_post_full_embedded
    );
    uintref4_type nu_post_prolong(
        solver.multigrid.nu_post_prolong.n = 0,
        solver.multigrid.nu_post_prolong.n_pre_full_embedded = 0,
        solver.multigrid.nu_post_prolong.n_full = 0,
        solver.multigrid.nu_post_prolong.n_post_full_embedded = 0
    );
    bool symmetrize_nu_pre_restrict_and_nu_post_prolong = true;

    std::string::const_iterator it = solver_str.begin();
    const bool b = qi::phrase_parse(
        it, solver_str.end(),
        qi::lit("petsc-cg")      [phx::ref(solver_id) = SOLVER_PARAMS::SOLVER_ID_PETSC_CG      ]
      | qi::lit("petsc-minres")  [phx::ref(solver_id) = SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES  ]
      | qi::lit("physbam-cg")    [phx::ref(solver_id) = SOLVER_PARAMS::SOLVER_ID_PHYSBAM_CG    ]
      | qi::lit("physbam-minres")[phx::ref(solver_id) = SOLVER_PARAMS::SOLVER_ID_PHYSBAM_MINRES]
      | (
            (qi::lit("mg(")   [phx::ref(solver_id) = SOLVER_PARAMS::SOLVER_ID_MG   ]
           | qi::lit("mgpcg(")[phx::ref(solver_id) = SOLVER_PARAMS::SOLVER_ID_MGPCG])
                >> -(smoother[phx::ref(smoother_id) = qi::_1] >> ',')
                >> -(qi::uint_[phx::ref(n_level) = qi::_1] >> ',')
                >> -(qi::uint_[phx::ref(mu) = qi::_1] >> ',')
                >> uint4[phx::ref(nu_pre_restrict) = qi::_1]
                >> -( ',' >>
                    uint4[
                        phx::ref(nu_post_prolong) = qi::_1,
                        phx::ref(symmetrize_nu_pre_restrict_and_nu_post_prolong) = false
                    ]
                )
            >> ')'
        ),
        qi::standard::blank
    );
    if(b) {
        if(it != solver_str.end())
            return 1;
        if(
            (solver_id == SOLVER_PARAMS::SOLVER_ID_MG ||
             solver_id == SOLVER_PARAMS::SOLVER_ID_MGPCG) &&
            symmetrize_nu_pre_restrict_and_nu_post_prolong
        )
            solver.multigrid.nu_post_prolong =
                solver.multigrid.nu_pre_restrict.Construct_Symmetric();
        return 0;
    }

    return 1;
}

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM
