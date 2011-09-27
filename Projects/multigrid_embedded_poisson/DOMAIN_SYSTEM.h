//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_SYSTEM_HPP

#include <cassert>

#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/container/vector/vector10.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/vector/vector10.hpp>
#include <boost/ref.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <Jeffrey_Utilities/ADD_REFERENCE_ADD_CONST.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>

#include "SYSTEM_SUM.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
class DOMAIN_SYSTEM
    : public SYSTEM_SUM< boost::mpl::vector2< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS > >
{
    typedef SYSTEM_SUM< boost::mpl::vector2< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS > > SYSTEM_SUM_;
public:
    typedef typename SYSTEM_SUM_::SCALAR_TYPE SCALAR_TYPE;
    static const int DIMENSION = SYSTEM_SUM_::DIMENSION;
    BOOST_MPL_ASSERT((boost::is_same<
        typename boost::remove_reference< T_REGULAR_SUBSYS >::type::SCALAR_TYPE,
        SCALAR_TYPE
    >));
    BOOST_MPL_ASSERT((boost::is_same<
        typename boost::remove_reference< T_EMBEDDING_SUBSYS >::type::SCALAR_TYPE,
        SCALAR_TYPE
    >));
    BOOST_MPL_ASSERT_RELATION( boost::remove_reference< T_REGULAR_SUBSYS >::type::DIMENSION, ==, DIMENSION );
    BOOST_MPL_ASSERT_RELATION( boost::remove_reference< T_EMBEDDING_SUBSYS >::type::DIMENSION, ==, DIMENSION );

    typedef typename SYSTEM_SUM_::MULTI_INDEX_TYPE MULTI_INDEX_TYPE;

    typedef T_REGULAR_SUBSYS REGULAR_SUBSYS_TYPE;
    typedef T_EMBEDDING_SUBSYS EMBEDDING_SUBSYS_TYPE;

    DOMAIN_SYSTEM(
        typename ADD_REFERENCE_ADD_CONST< T_REGULAR_SUBSYS >::type regular_subsys,
        typename ADD_REFERENCE_ADD_CONST< T_EMBEDDING_SUBSYS >::type embedding_subsys);

    typename boost::add_reference< T_REGULAR_SUBSYS >::type Regular_Subsys();
    typename ADD_REFERENCE_ADD_CONST< T_REGULAR_SUBSYS >::type Regular_Subsys() const;
    typename boost::add_reference< T_EMBEDDING_SUBSYS >::type Embedding_Subsys();
    typename ADD_REFERENCE_ADD_CONST< T_EMBEDDING_SUBSYS >::type Embedding_Subsys() const;

    typedef MULTI_INDEX_BOUND< DIMENSION > MULTI_INDEX_BOUND_TYPE;
    MULTI_INDEX_BOUND_TYPE Multi_Index_Bound() const;
    MULTI_INDEX_BOUND_TYPE Cell_Multi_Index_Bound() const;
    typedef VECTOR< SCALAR_TYPE, DIMENSION > VECTOR_TYPE;
    VECTOR_TYPE Dx() const;

    template< class T_RHS_OF_INDEX >
    void Set_Dirichlet_Grid_BC(
        const MULTI_INDEX_TYPE& multi_index,
        const SCALAR_TYPE p,
        const T_RHS_OF_INDEX& rhs_of_index);

    using SYSTEM_SUM_::Diag;
    using SYSTEM_SUM_::Stencil_N_Nonzero;
    using SYSTEM_SUM_::Stencil_Sum;
    using SYSTEM_SUM_::Set_Dirichlet_Grid_BC;
    using SYSTEM_SUM_::Apply;
    using SYSTEM_SUM_::Apply_Transpose;
};

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
inline DOMAIN_SYSTEM< T_REGULAR_SUBSYS&, T_EMBEDDING_SUBSYS& >
Make_Domain_System(T_REGULAR_SUBSYS& regular_subsys, T_EMBEDDING_SUBSYS& embedding_subsys)
{ return DOMAIN_SYSTEM< T_REGULAR_SUBSYS&, T_EMBEDDING_SUBSYS& >(regular_subsys, embedding_subsys); }

//#####################################################################
//#####################################################################

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
inline 
DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::
DOMAIN_SYSTEM(
    typename ADD_REFERENCE_ADD_CONST< T_REGULAR_SUBSYS >::type regular_subsys,
    typename ADD_REFERENCE_ADD_CONST< T_EMBEDDING_SUBSYS >::type embedding_subsys)
    : SYSTEM_SUM_(boost::fusion::make_vector(boost::ref(regular_subsys), boost::ref(embedding_subsys)))
{ }

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
inline typename boost::add_reference< T_REGULAR_SUBSYS >::type
DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::
Regular_Subsys()
{ return boost::fusion::at_c<0>(SYSTEM_SUM_::systems); }

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
inline typename ADD_REFERENCE_ADD_CONST< T_REGULAR_SUBSYS >::type
DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::
Regular_Subsys() const
{ return boost::fusion::at_c<0>(SYSTEM_SUM_::systems); }

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
inline typename boost::add_reference< T_EMBEDDING_SUBSYS >::type
DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::
Embedding_Subsys()
{ return boost::fusion::at_c<1>(SYSTEM_SUM_::systems); }

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
inline typename ADD_REFERENCE_ADD_CONST< T_EMBEDDING_SUBSYS >::type
DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::
Embedding_Subsys() const
{ return boost::fusion::at_c<1>(SYSTEM_SUM_::systems); }

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
inline typename DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::MULTI_INDEX_BOUND_TYPE
DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::
Multi_Index_Bound() const
{
    assert(Regular_Subsys().multi_index_bound == Embedding_Subsys().multi_index_bound);
    return Regular_Subsys().multi_index_bound;
}

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
inline typename DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::MULTI_INDEX_BOUND_TYPE
DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::
Cell_Multi_Index_Bound() const
{ return Regular_Subsys().Cell_Multi_Index_Bound(); }

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
inline typename DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::VECTOR_TYPE
DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::
Dx() const
{ return Regular_Subsys().dx; }

template< class T_REGULAR_SUBSYS, class T_EMBEDDING_SUBSYS >
template< class T_RHS_OF_INDEX >
inline void
DOMAIN_SYSTEM< T_REGULAR_SUBSYS, T_EMBEDDING_SUBSYS >::
Set_Dirichlet_Grid_BC(
    const MULTI_INDEX_TYPE& multi_index,
    const SCALAR_TYPE p,
    const T_RHS_OF_INDEX& rhs_of_index)
{
    const int linear_index = Multi_Index_Bound().Linear_Index(multi_index);
    Set_Dirichlet_Grid_BC(linear_index, p, rhs_of_index);
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_DOMAIN_SYSTEM_HPP
