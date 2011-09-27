//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SYSTEM_SUM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SYSTEM_SUM_HPP

#include <cassert>

#include <functional>

#include <boost/concept/assert.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/utility/result_of.hpp>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T_SYSTEM_SEQUENCE >
struct SYSTEM_SUM
{
    BOOST_MPL_ASSERT((boost::mpl::is_sequence< T_SYSTEM_SEQUENCE >));
    typedef typename boost::remove_reference<
        typename boost::mpl::front< T_SYSTEM_SEQUENCE >::type
    >::type::SCALAR_TYPE SCALAR_TYPE;
    static const int DIMENSION = boost::remove_reference<
        typename boost::mpl::front< T_SYSTEM_SEQUENCE >::type
    >::type::DIMENSION;

    typedef VECTOR< int, DIMENSION > MULTI_INDEX_TYPE;

    typedef T_SYSTEM_SEQUENCE SYSTEM_SEQUENCE_TYPE;

    typename boost::fusion::result_of::as_vector< T_SYSTEM_SEQUENCE >::type systems;

    template< class T_SYSTEM_SEQUENCE2 >
    explicit SYSTEM_SUM(const T_SYSTEM_SEQUENCE2& systems_,
        typename boost::enable_if< boost::is_convertible<
            const T_SYSTEM_SEQUENCE2&,
            typename boost::fusion::result_of::as_vector< T_SYSTEM_SEQUENCE >::type
        > >::type* = 0);

    SCALAR_TYPE Diag(const int index) const;
    // This is really just an upper bound...
    int Stencil_N_Nonzero(const int index) const;
    SCALAR_TYPE Stencil_Sum(const int index) const;

    template< class T_STENCIL_PROXY >
    void Add_Stencil_To(
        const int index,
        const T_STENCIL_PROXY& stencil_proxy) const;

    template< class T_RHS_OF_INDEX >
    void Set_Dirichlet_Grid_BC(
        const int index,
        const SCALAR_TYPE p,
        const T_RHS_OF_INDEX& rhs_of_index);

    void Apply(const ARRAY_VIEW< const SCALAR_TYPE > x, ARRAY_VIEW< SCALAR_TYPE > y) const;
    void Apply(const ARRAY_VIEW< const SCALAR_TYPE > x, ARRAY_VIEW< SCALAR_TYPE > y, const MULTI_INDEX_TYPE& strides) const;
    SCALAR_TYPE Apply(const int index, const ARRAY_VIEW< const SCALAR_TYPE > x) const;
    SCALAR_TYPE Apply(const int index, const ARRAY_VIEW< const SCALAR_TYPE > x, const MULTI_INDEX_TYPE& strides) const;
    void Apply_Transpose(const int index, ARRAY_VIEW< SCALAR_TYPE > y, const SCALAR_TYPE x) const;
    void Apply_Transpose(const int index, ARRAY_VIEW< SCALAR_TYPE > y, const SCALAR_TYPE x, const MULTI_INDEX_TYPE& strides) const;
};

//#####################################################################
//#####################################################################

template< class T_SYSTEM_SEQUENCE >
template< class T_SYSTEM_SEQUENCE2 >
inline
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
SYSTEM_SUM(const T_SYSTEM_SEQUENCE2& systems_,
    typename boost::enable_if< boost::is_convertible<
        const T_SYSTEM_SEQUENCE2&,
        typename boost::fusion::result_of::as_vector< T_SYSTEM_SEQUENCE >::type
    > >::type*)
    : systems(systems_)
{ }

namespace Detail_SYSTEM_SUM
{

template< class T >
struct DIAG
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        DIAG, (( int const, index ))
    )
public:
    typedef T result_type;
    template< class T_SYSTEM >
    T operator()(const T_SYSTEM& system) const
    { return system.Valid_Index(index) ? system.Diag(index) : static_cast<T>(0); }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
inline typename SYSTEM_SUM< T_SYSTEM_SEQUENCE >::SCALAR_TYPE
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Diag(const int index) const
{
    typedef Detail_SYSTEM_SUM::DIAG< SCALAR_TYPE > DIAG_;
    return boost::fusion::accumulate(
        boost::fusion::transform(systems, DIAG_(index)),
        static_cast< SCALAR_TYPE >(0),
        std::plus< SCALAR_TYPE >()
    );
}

namespace Detail_SYSTEM_SUM
{

struct STENCIL_N_NONZERO
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        STENCIL_N_NONZERO, (( int const, index ))
    )
public:
    typedef int result_type;
    template< class T_SYSTEM >
    int operator()(const T_SYSTEM& system) const
    { return system.Valid_Index(index) ? system.Stencil_N_Nonzero(index) : 0; }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
inline int
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Stencil_N_Nonzero(const int index) const
{
    typedef Detail_SYSTEM_SUM::STENCIL_N_NONZERO STENCIL_N_NONZERO_;
    return boost::fusion::accumulate(
        boost::fusion::transform(systems, STENCIL_N_NONZERO_(index)),
        0,
        std::plus< int >()
    );
}

namespace Detail_SYSTEM_SUM
{

template< class T >
struct STENCIL_SUM
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        STENCIL_SUM, (( int const, index ))
    )
public:
    typedef T result_type;
    template< class T_SYSTEM >
    T operator()(const T_SYSTEM& system) const
    { return system.Valid_Index(index) ? system.Stencil_Sum(index) : static_cast<T>(0); }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
inline typename SYSTEM_SUM< T_SYSTEM_SEQUENCE >::SCALAR_TYPE
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Stencil_Sum(const int index) const
{
    typedef Detail_SYSTEM_SUM::STENCIL_SUM< SCALAR_TYPE > STENCIL_SUM_;
    return boost::fusion::accumulate(
        boost::fusion::transform(systems, STENCIL_SUM_(index)),
        static_cast< SCALAR_TYPE >(0),
        std::plus< SCALAR_TYPE >()
    );
}

namespace Detail_SYSTEM_SUM
{

template< class T_STENCIL_PROXY >
struct ADD_STENCIL_TO
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        ADD_STENCIL_TO,
        (( /******/ int const, index ))
        (( typename T_STENCIL_PROXY const &, stencil_proxy ))
    )
public:
    typedef void result_type;
    template< class T_SYSTEM >
    void operator()(const T_SYSTEM& system) const
    {
        if(system.Valid_Index(index))
            system.Add_Stencil_To(index, stencil_proxy);
    }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
template< class T_STENCIL_PROXY >
inline void
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Add_Stencil_To(
    const int index,
    const T_STENCIL_PROXY& stencil_proxy) const
{
    BOOST_CONCEPT_ASSERT((STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    typedef Detail_SYSTEM_SUM::ADD_STENCIL_TO< T_STENCIL_PROXY > ADD_STENCIL_TO_;
    boost::fusion::for_each(systems, ADD_STENCIL_TO_(index, stencil_proxy));
}

namespace Detail_SYSTEM_SUM
{

template< class T, class T_RHS_OF_INDEX >
struct SET_DIRICHLET_GRID_BC
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SET_DIRICHLET_GRID_BC,
        (( /******/ int const, index ))
        (( typename T const, p ))
        (( typename T_RHS_OF_INDEX const, rhs_of_index ))
    )
public:
    typedef void result_type;
    template< class T_SYSTEM >
    void operator()(T_SYSTEM& system) const
    {
        if(system.Valid_Index(index))
            system.Set_Dirichlet_Grid_BC(index, p, rhs_of_index);
    }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
template< class T_RHS_OF_INDEX >
inline void
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Set_Dirichlet_Grid_BC(
    const int index,
    const SCALAR_TYPE p,
    const T_RHS_OF_INDEX& rhs_of_index)
{
    typedef Detail_SYSTEM_SUM::SET_DIRICHLET_GRID_BC< SCALAR_TYPE, T_RHS_OF_INDEX > SET_DIRICHLET_GRID_BC_;
    boost::fusion::for_each(systems, SET_DIRICHLET_GRID_BC_(index, p, rhs_of_index));
}

namespace Detail_SYSTEM_SUM
{

template< class T >
struct APPLY
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        APPLY,
        (( typename ARRAY_VIEW<const T> const, x ))
        (( typename ARRAY_VIEW<T>&, y ))
    )
public:
    typedef void result_type;
    template< class T_SYSTEM >
    void operator()(const T_SYSTEM& system) const
    { system.Apply(x, y); }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
inline void
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Apply(const ARRAY_VIEW< const SCALAR_TYPE > x, ARRAY_VIEW< SCALAR_TYPE > y) const
{
    typedef Detail_SYSTEM_SUM::APPLY< SCALAR_TYPE > APPLY_;
    boost::fusion::for_each(systems, APPLY_(x, y));
}

namespace Detail_SYSTEM_SUM
{

template< class T, int D >
struct APPLY_WITH_STRIDES
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        APPLY_WITH_STRIDES,
        (( typename ARRAY_VIEW<const T> const, x ))
        (( typename ARRAY_VIEW<T>&, y ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<int,D> )) const, strides ))
    )
public:
    typedef void result_type;
    template< class T_SYSTEM >
    void operator()(const T_SYSTEM& system) const
    { system.Apply(x, y, strides); }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
inline void
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Apply(const ARRAY_VIEW< const SCALAR_TYPE > x, ARRAY_VIEW< SCALAR_TYPE > y, const MULTI_INDEX_TYPE& strides) const
{
    typedef Detail_SYSTEM_SUM::APPLY_WITH_STRIDES< SCALAR_TYPE, DIMENSION > APPLY_WITH_STRIDES_;
    boost::fusion::for_each(systems, APPLY_WITH_STRIDES_(x, y, strides));
}

namespace Detail_SYSTEM_SUM
{

template< class T >
struct APPLY_AT
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        APPLY_AT,
        (( /******/ int const, index ))
        (( typename ARRAY_VIEW<const T> const, x ))
    )
public:
    typedef T result_type;
    template< class T_SYSTEM >
    T operator()(const T_SYSTEM& system) const
    { return system.Valid_Index(index) ? system.Apply(index, x) : static_cast<T>(0); }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
inline typename SYSTEM_SUM< T_SYSTEM_SEQUENCE >::SCALAR_TYPE
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Apply(const int index, const ARRAY_VIEW< const SCALAR_TYPE > x) const
{
    typedef Detail_SYSTEM_SUM::APPLY_AT< SCALAR_TYPE > APPLY_AT_;
    return boost::fusion::accumulate(
        boost::fusion::transform(systems, APPLY_AT_(index, x)),
        static_cast< SCALAR_TYPE >(0),
        std::plus< SCALAR_TYPE >()
    );
}

namespace Detail_SYSTEM_SUM
{

template< class T, int D >
struct APPLY_AT_WITH_STRIDES
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        APPLY_AT_WITH_STRIDES,
        (( /******/ int const, index ))
        (( typename ARRAY_VIEW<const T> const, x ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<int,D> )) const, strides ))
    )
public:
    typedef T result_type;
    template< class T_SYSTEM >
    T operator()(const T_SYSTEM& system) const
    { return system.Valid_Index(index) ? system.Apply(index, x, strides) : static_cast<T>(0); }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
inline typename SYSTEM_SUM< T_SYSTEM_SEQUENCE >::SCALAR_TYPE
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Apply(const int index, const ARRAY_VIEW< const SCALAR_TYPE > x, const MULTI_INDEX_TYPE& strides) const
{
    typedef Detail_SYSTEM_SUM::APPLY_AT_WITH_STRIDES< SCALAR_TYPE, DIMENSION > APPLY_AT_WITH_STRIDES_;
    return boost::fusion::accumulate(
        boost::fusion::transform(systems, APPLY_AT_WITH_STRIDES_(index, x, strides)),
        static_cast< SCALAR_TYPE >(0),
        std::plus< SCALAR_TYPE >()
    );
}

namespace Detail_SYSTEM_SUM
{

template< class T >
struct APPLY_TRANSPOSE
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        APPLY_TRANSPOSE,
        (( /******/ int const, index ))
        (( typename ARRAY_VIEW<T>&, y ))
        (( typename T const, x ))
    )
public:
    typedef void result_type;
    template< class T_SYSTEM >
    void operator()(const T_SYSTEM& system) const
    {
        if(system.Valid_Index(index))
            system.Apply_Transpose(index, y, x);
    }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
inline void
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Apply_Transpose(const int index, ARRAY_VIEW< SCALAR_TYPE > y, const SCALAR_TYPE x) const
{
    typedef Detail_SYSTEM_SUM::APPLY_TRANSPOSE< SCALAR_TYPE > APPLY_TRANSPOSE_;
    boost::fusion::for_each(systems, APPLY_TRANSPOSE_(index, y, x));
}

namespace Detail_SYSTEM_SUM
{

template< class T, int D >
struct APPLY_TRANSPOSE_WITH_STRIDES
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        APPLY_TRANSPOSE_WITH_STRIDES,
        (( /******/ int const, index ))
        (( typename ARRAY_VIEW<T>&, y ))
        (( typename T const, x ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( VECTOR<int,D> )) const, strides ))
    )
public:
    typedef void result_type;
    template< class T_SYSTEM >
    void operator()(const T_SYSTEM& system) const
    {
        if(system.Valid_Index(index))
            system.Apply_Tranpose(index, y, x, strides);
    }
};

} // namespace Detail_SYSTEM_SUM

template< class T_SYSTEM_SEQUENCE >
inline void
SYSTEM_SUM< T_SYSTEM_SEQUENCE >::
Apply_Transpose(const int index, ARRAY_VIEW< SCALAR_TYPE > y, const SCALAR_TYPE x, const MULTI_INDEX_TYPE& strides) const
{
    typedef Detail_SYSTEM_SUM::APPLY_TRANSPOSE_WITH_STRIDES< SCALAR_TYPE, DIMENSION > APPLY_TRANSPOSE_WITH_STRIDES_;
    boost::fusion::for_each(systems, APPLY_TRANSPOSE_WITH_STRIDES_(index, y, x, strides));
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_SYSTEM_SUM_HPP
