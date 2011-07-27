//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CUBE2_SIMPLEX_PARTITION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CUBE2_SIMPLEX_PARTITION_HPP

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/Geometry/CUBE2_SIMPLEX_PARTITION_ITERATOR.h>

namespace PhysBAM
{

template< int D >
struct CUBE2_SIMPLEX_PARTITION
{
    static const int CENTER_INDEX = CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::CENTER_INDEX;

    typedef CUBE2_SIMPLEX_PARTITION_ITERATOR<D> iterator;
    typedef iterator const_iterator;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

template< int D >
inline typename CUBE2_SIMPLEX_PARTITION<D>::iterator
CUBE2_SIMPLEX_PARTITION<D>::
begin() const
{ return iterator(BEGIN_TAG()); }

template< int D >
inline typename CUBE2_SIMPLEX_PARTITION<D>::iterator
CUBE2_SIMPLEX_PARTITION<D>::
end() const
{ return iterator(END_TAG()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CUBE2_SIMPLEX_PARTITION_HPP
