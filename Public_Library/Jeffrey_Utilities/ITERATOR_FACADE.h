//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ITERATOR_FACADE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ITERATOR_FACADE_HPP

#include <cstddef>

#include <boost/iterator/iterator_facade.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

namespace Detail_ITERATOR_FACADE
{

template< class T_REFERENCE >
struct OPERATOR_ARROW_DISPATCH;

} // namespace Detail_ITERATOR_FACADE

template<
    class Derived,
    class Value,
    class CategoryOrTraversal,
    class Reference = Value&,
    class Difference = std::ptrdiff_t
>
class ITERATOR_FACADE
    : public boost::iterator_facade<
          Derived,
          Value,
          CategoryOrTraversal,
          Reference,
          Difference
      >
{
    typedef boost::iterator_facade<
        Derived,
        Value,
        CategoryOrTraversal,
        Reference,
        Difference
    > iterator_facade_;
    typedef Detail_ITERATOR_FACADE::OPERATOR_ARROW_DISPATCH< Reference > OPERATOR_ARROW_DISPATCH_;
public:
    typedef typename OPERATOR_ARROW_DISPATCH_::result_type pointer;
    pointer operator->() const
    { return OPERATOR_ARROW_DISPATCH_::Apply(iterator_facade_::operator*()); }
};

namespace Detail_ITERATOR_FACADE
{

template< class T_REFERENCE >
struct OPERATOR_ARROW_DISPATCH
{
    class PROXY
    {
        PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
            PROXY, (( typename T_REFERENCE const, ref ))
        )
        friend struct OPERATOR_ARROW_DISPATCH;
    public:
        const T_REFERENCE* operator->() const
        { return &ref; }
    };
    typedef PROXY result_type;
    static result_type Apply(const T_REFERENCE& x)
    { return result_type(x); }
};

template< class T >
struct OPERATOR_ARROW_DISPATCH< T& >
{
    typedef T* result_type;
    static result_type Apply(T& x)
    { return &x; }
};

} // namespace Detail_ITERATOR_FACADE

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ITERATOR_FACADE_HPP
