//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ITERATOR_ADAPTOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ITERATOR_ADAPTOR_HPP

#include <boost/iterator/iterator_adaptor.hpp>

#include <Jeffrey_Utilities/ITERATOR_FACADE.h>

namespace PhysBAM
{

template<
    class Derived,
    class Base,
    class Value               = boost::use_default,
    class CategoryOrTraversal = boost::use_default,
    class Reference           = boost::use_default,
    class Difference          = boost::use_default
>
class ITERATOR_ADAPTOR
    : public boost::iterator_adaptor<
          Derived,
          Base,
          Value,
          CategoryOrTraversal,
          Reference,
          Difference
      >
{
    typedef boost::iterator_adaptor<
        Derived,
        Base,
        Value,
        CategoryOrTraversal,
        Reference,
        Difference
    > iterator_adaptor_;
    typedef Detail_ITERATOR_FACADE::OPERATOR_ARROW_DISPATCH<
        typename iterator_adaptor_::reference
    > OPERATOR_ARROW_DISPATCH_;
public:

    ITERATOR_ADAPTOR()
    { }

    explicit ITERATOR_ADAPTOR(const Base& base)
        : iterator_adaptor_(base)
    { }

    typedef typename OPERATOR_ARROW_DISPATCH_::result_type pointer;
    pointer operator->() const
    { return OPERATOR_ARROW_DISPATCH_::Apply(iterator_adaptor_::operator*()); }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ITERATOR_ADAPTOR_HPP
