//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <string>

#include <boost/function.hpp>

#include <boost/fusion/container/vector/vector10.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/spirit/home/phoenix/core/reference.hpp>
#include <boost/spirit/home/phoenix/fusion/at.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/spirit/home/phoenix/statement/sequence.hpp>
#include <boost/spirit/home/qi/action/action.hpp>
#include <boost/spirit/home/qi/auxiliary/eps.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/char/char_class.hpp>
#include <boost/spirit/home/qi/nonterminal/rule.hpp>
#include <boost/spirit/home/qi/numeric/real.hpp>
#include <boost/spirit/home/qi/numeric/uint.hpp>
#include <boost/spirit/home/qi/operator/alternative.hpp>
#include <boost/spirit/home/qi/operator/optional.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/parse.hpp>
#include <boost/spirit/home/qi/parse_attr.hpp>
#include <boost/spirit/home/qi/string/lit.hpp>
#include <boost/spirit/home/support/argument.hpp>
#include <boost/spirit/home/support/context.hpp>

#include <Jeffrey_Utilities/Functional/Post_Compose_Scale.h>
#include <Jeffrey_Utilities/Functional/Post_Compose_Translate.h>

#include "../Main_Examples.h"
#include "../Params/EXAMPLE_PARAMS.h"

#include "Parse_String_Example.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace
{

template< class T, int D >
inline void Get_Main_Example_Helper(
    const boost::fusion::vector3< unsigned int, float, float > u_id_transform,
    const boost::fusion::vector2< unsigned int, float > beta_id_transform,
    boost::function< T ( const VECTOR<T,D>& ) >& u,
    boost::function< VECTOR<T,D> ( const VECTOR<T,D>& ) >& grad_u,
    boost::function< T ( const VECTOR<T,D>& ) >& beta,
    boost::function< T ( const VECTOR<T,D>& ) >& f);

} // namespace

template< class T, int D >
int Parse_String_Example(
    const std::string& example_str,
    typename EXAMPLE_PARAMS<T,D>::PROBLEM_TYPE& problem)
{
    namespace phx = boost::phoenix;
    namespace qi  = boost::spirit::qi;

    typedef EXAMPLE_PARAMS_BASE::CONSTRAINT_ID CONSTRAINT_ID;

    static const qi::rule<
        std::string::const_iterator,
        boost::fusion::vector3< unsigned int, float, float > ( ),
        qi::standard::blank_type
    > u_id_transform_ =
        qi::eps [
            phx::at_c<1>(qi::_val) = static_cast< float >(1),
            phx::at_c<2>(qi::_val) = static_cast< float >(0)
        ]
     >> -(qi::float_ >> '*')       [ phx::at_c<1>(qi::_val) = qi::_1 ]
     >>  ('{' >> qi::uint_ >> '}') [ phx::at_c<0>(qi::_val) = qi::_1 ]
     >> -('+' >> qi::float_)       [ phx::at_c<2>(qi::_val) = qi::_1 ];
 
    static const qi::rule<
        std::string::const_iterator,
        boost::fusion::vector2< unsigned int, float > ( ),
        qi::standard::blank_type
    > beta_id_transform_ =
        qi::eps [ phx::at_c<1>(qi::_val) = static_cast< float >(1) ]
     >> -(qi::float_ >> '*')      [ phx::at_c<1>(qi::_val) = qi::_1 ]
     >> ('{' >> qi::uint_ >> '}') [ phx::at_c<0>(qi::_val) = qi::_1 ];

    static const qi::rule<
        std::string::const_iterator,
        boost::fusion::vector2< CONSTRAINT_ID, float > ( ),
        qi::standard::blank_type
    > constraint_ =
        qi::eps [ phx::at_c<1>(qi::_val) = static_cast< float >(0) ]
     >> qi::lit("single-cell") [ phx::at_c<0>(qi::_val) = EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_SINGLE_CELL ]
      | (
            (qi::lit("double-cell") [ phx::at_c<0>(qi::_val) = EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_DOUBLE_CELL ]
           | qi::lit("aggregate")   [ phx::at_c<0>(qi::_val) = EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_AGGREGATE   ])
            >> -( '(' >> qi::float_ >> ')' ) [ phx::at_c<1>(qi::_val) = qi::_1 ]
        );

    std::string::const_iterator it;

    {
        boost::fusion::vector3< unsigned int, float, float > u_id_transform;
        boost::fusion::vector2< unsigned int, float > beta_id_transform;
        const bool b = qi::phrase_parse(
            it = example_str.begin(), example_str.end(),
            qi::lit("neumann(") >> u_id_transform_ >> ',' >> beta_id_transform_ >> ')',
            qi::standard::blank,
            u_id_transform, beta_id_transform
        );
        if(b) {
            if(it != example_str.end())
                return 1;
            typename EXAMPLE_PARAMS<T,D>::NEUMANN_PARAMS neumann_problem;
            Get_Main_Example_Helper<T,D>(
                u_id_transform, beta_id_transform,
                neumann_problem.u,    neumann_problem.grad_u,
                neumann_problem.beta, neumann_problem.f
            );
            problem = neumann_problem;
            return 0;
        }
    }

    {
        boost::fusion::vector3< unsigned int, float, float > u_id_transform;
        boost::fusion::vector2< unsigned int, float > beta_id_transform;
        CONSTRAINT_ID constraint_id;
        float min_relative_indy_weight = 0;
        const bool b = qi::phrase_parse(
            it = example_str.begin(), example_str.end(),
            qi::lit("dirichlet(")
                >>    u_id_transform_[phx::ref(   u_id_transform) = qi::_1] >> ','
                >> beta_id_transform_[phx::ref(beta_id_transform) = qi::_1] >> ','
                >> constraint_[
                    phx::ref(constraint_id) = phx::at_c<0>(qi::_1),
                    phx::ref(min_relative_indy_weight) = phx::at_c<1>(qi::_1)
                ]
            >> ')',
            qi::standard::blank
        );
        if(b) {
            if(it != example_str.end())
                return 1;
            typename EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS dirichlet_problem;
            dirichlet_problem.constraint_id = constraint_id;
            dirichlet_problem.min_relative_indy_weight = min_relative_indy_weight;
            Get_Main_Example_Helper<T,D>(
                u_id_transform, beta_id_transform,
                dirichlet_problem.u,    dirichlet_problem.grad_u,
                dirichlet_problem.beta, dirichlet_problem.f
            );
            problem = dirichlet_problem;
            return 0;
        }
    }

    {
        boost::fusion::vector3< unsigned int, float, float > u_m_id_transform;
        boost::fusion::vector2< unsigned int, float > beta_m_id_transform;
        boost::fusion::vector3< unsigned int, float, float > u_p_id_transform;
        boost::fusion::vector2< unsigned int, float > beta_p_id_transform;
        CONSTRAINT_ID constraint_id;
        float min_relative_indy_weight = 0;
        const bool b = qi::phrase_parse(
            it = example_str.begin(), example_str.end(),
            qi::lit("interface(")
                >> '('
                    >>    u_id_transform_[phx::ref(   u_m_id_transform) = qi::_1] >> ','
                    >> beta_id_transform_[phx::ref(beta_m_id_transform) = qi::_1]
                >> ')'
                >> ','
                >> '('
                    >>    u_id_transform_[phx::ref(   u_p_id_transform) = qi::_1] >> ','
                    >> beta_id_transform_[phx::ref(beta_p_id_transform) = qi::_1]
                >> ')'
                >> ','
                >> constraint_[
                    phx::ref(constraint_id) = phx::at_c<0>(qi::_1),
                    phx::ref(min_relative_indy_weight) = phx::at_c<1>(qi::_1)
                ]
            >> ')',
            qi::standard::blank
        );
        if(b) {
            if(it != example_str.end())
                return 1;
            typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS interface_problem;
            interface_problem.constraint_id = constraint_id;
            interface_problem.min_relative_indy_weight = min_relative_indy_weight;
            Get_Main_Example_Helper<T,D>(
                u_m_id_transform, beta_m_id_transform,
                interface_problem.negative.u,    interface_problem.negative.grad_u,
                interface_problem.negative.beta, interface_problem.negative.f
            );
            Get_Main_Example_Helper<T,D>(
                u_p_id_transform, beta_p_id_transform,
                interface_problem.positive.u,    interface_problem.positive.grad_u,
                interface_problem.positive.beta, interface_problem.positive.f
            );
            problem = interface_problem;
            return 0;
        }
    }

    return 1;
}

namespace
{

template< class T, int D >
inline void
Get_Main_Example_Helper(
    const boost::fusion::vector3< unsigned int, float, float > u_id_transform,
    const boost::fusion::vector2< unsigned int, float > beta_id_transform,
    boost::function< T ( const VECTOR<T,D>& ) >& u,
    boost::function< VECTOR<T,D> ( const VECTOR<T,D>& ) >& grad_u,
    boost::function< T ( const VECTOR<T,D>& ) >& beta,
    boost::function< T ( const VECTOR<T,D>& ) >& f)
{
    const unsigned int    u_id        = boost::fusion::at_c<0>(u_id_transform);
    const float           u_scale     = boost::fusion::at_c<1>(u_id_transform);
    const float           u_translate = boost::fusion::at_c<2>(u_id_transform);
    const unsigned int beta_id        = boost::fusion::at_c<0>(beta_id_transform);
    const float        beta_scale     = boost::fusion::at_c<1>(beta_id_transform);
    Get_Main_Example<T,D>(u_id, beta_id, u, grad_u, beta, f);
    const unsigned int case_ =
        ((   u_scale     != 1) << 0)
      | ((   u_translate != 0) << 1)
      | ((beta_scale     != 1) << 2);
    switch(case_) {
    case 0:
        break;
    case 1:
        u = Post_Compose_Scale(u_scale, u);
        f = Post_Compose_Scale(u_scale, f);
        break;
    case 2:
        u = Post_Compose_Translate(u_translate, u);
        break;
    case 3:
        u = Post_Compose_Translate(
            u_translate,
            Post_Compose_Scale(u_scale, u)
        );
        f = Post_Compose_Scale(u_scale, f);
        break;
    case 4:
        beta = Post_Compose_Scale(beta_scale, beta);
        f = Post_Compose_Scale(beta_scale, beta);
        break;
    case 5:
        u = Post_Compose_Scale(u_scale, u);
        beta = Post_Compose_Scale(beta_scale, beta);
        f = Post_Compose_Scale(u_scale * beta_scale, f);
        break;
    case 6:
        u = Post_Compose_Translate(u_translate, u);
        beta = Post_Compose_Scale(beta_scale, beta);
        f = Post_Compose_Scale(beta_scale, beta);
        break;
    case 7:
        u = Post_Compose_Translate(
            u_translate,
            Post_Compose_Scale(u_scale, u)
        );
        beta = Post_Compose_Scale(beta_scale, beta);
        f = Post_Compose_Scale(u_scale * beta_scale, f);
        break;
    }
}

} // namespace

#define EXPLICIT_INSTANTIATION( T, D ) \
template int Parse_String_Example<T,D>( \
    const std::string& example_str, \
    EXAMPLE_PARAMS<T,D>::PROBLEM_TYPE& problem);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
