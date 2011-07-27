//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARAMS_EXAMPLE_PARAMS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARAMS_EXAMPLE_PARAMS_HPP

#include <boost/blank.hpp>
#include <boost/function.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/variant.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

struct EXAMPLE_PARAMS_BASE
{
    enum GRID_BC_ID
    {
        GRID_BC_ID_NULL = 0,
        GRID_BC_ID_NEUMANN,
        GRID_BC_ID_DIRICHLET
    } grid_bc_id;

    enum CONSTRAINT_ID
    {
        CONSTRAINT_ID_NULL = 0,
        CONSTRAINT_ID_SINGLE_CELL,
        CONSTRAINT_ID_DOUBLE_CELL,
        CONSTRAINT_ID_AGGREGATE
    };

    EXAMPLE_PARAMS_BASE()
        : grid_bc_id(GRID_BC_ID_NULL)
    { }
};

template< class T, int D >
struct EXAMPLE_PARAMS
    : EXAMPLE_PARAMS_BASE
{
    typedef EXAMPLE_PARAMS_BASE::CONSTRAINT_ID CONSTRAINT_ID;

    //typedef T (*SCALAR_FUNCTION_TYPE)( const VECTOR<T,D>& );
    //typedef VECTOR<T,D> (*VECTOR_FUNCTION_TYPE)( const VECTOR<T,D>& );
    typedef boost::function< T ( const VECTOR<T,D>& ) > SCALAR_FUNCTION_TYPE;
    typedef boost::function< VECTOR<T,D> ( const VECTOR<T,D>& ) > VECTOR_FUNCTION_TYPE;

    struct DOMAIN_PARAMS
    {
        SCALAR_FUNCTION_TYPE u;
        VECTOR_FUNCTION_TYPE grad_u;
        SCALAR_FUNCTION_TYPE beta;
        SCALAR_FUNCTION_TYPE f;
        //DOMAIN_PARAMS()
        //    : u(0), grad_u(0), beta(0), f(0)
        //{ }
    };

    struct NEUMANN_PARAMS
        : DOMAIN_PARAMS
    { };

    struct DIRICHLET_PARAMS
        : DOMAIN_PARAMS
    {
        typedef EXAMPLE_PARAMS_BASE::CONSTRAINT_ID CONSTRAINT_ID;
        CONSTRAINT_ID constraint_id;
        float min_relative_indy_weight;
        DIRICHLET_PARAMS()
            : constraint_id(CONSTRAINT_ID_NULL),
              min_relative_indy_weight(0)
        { }
    };

    struct INTERFACE_PARAMS
    {
        typedef EXAMPLE_PARAMS_BASE::CONSTRAINT_ID CONSTRAINT_ID;
        CONSTRAINT_ID constraint_id;
        float min_relative_indy_weight;
        DOMAIN_PARAMS negative;
        DOMAIN_PARAMS positive;
        INTERFACE_PARAMS()
            : constraint_id(CONSTRAINT_ID_NULL),
              min_relative_indy_weight(0)
        { }
    };

    typedef boost::variant<
        boost::blank,
        NEUMANN_PARAMS, DIRICHLET_PARAMS, INTERFACE_PARAMS
    > PROBLEM_TYPE;
    PROBLEM_TYPE problem;

    bool Verify() const;

private:
    struct VERIFY_VISITOR;
};

//#####################################################################
//#####################################################################

template< class T, int D >
struct EXAMPLE_PARAMS<T,D>::
VERIFY_VISITOR
{
    typedef bool result_type;
    bool operator()(boost::blank) const
    { return false; }
    bool operator()(const DOMAIN_PARAMS& problem) const
    { return problem.u && problem.grad_u && problem.beta && problem.f; }
    bool operator()(const INTERFACE_PARAMS& problem) const
    { return operator()(problem.negative) && operator()(problem.positive); }
};

template< class T, int D >
inline bool
EXAMPLE_PARAMS<T,D>::
Verify() const
{ return boost::apply_visitor(VERIFY_VISITOR(), problem); }

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARAMS_EXAMPLE_PARAMS_HPP
