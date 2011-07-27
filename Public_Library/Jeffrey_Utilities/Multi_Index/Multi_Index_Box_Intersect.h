//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOX_INTERSECT_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOX_INTERSECT_HPP

#include <boost/mpl/assert.hpp>

//#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOX.h>
//#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
//#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>

namespace PhysBAM
{

template< class BOX1, class BOX2 >
struct MULTI_INDEX_BOX_INTERSECT_RESULT
{
    BOOST_MPL_ASSERT_RELATION( BOX1::DIMENSION, ==, BOX2::DIMENSION );
    typedef MULTI_INDEX_BOX< BOX1::DIMENSION > type;
};

template< class BOX1, class BOX2 >
inline MULTI_INDEX_BOX< BOX1::DIMENSION >
Multi_Index_Box_Intersect(const BOX1& box1, const BOX2& box2)
{
    BOOST_MPL_ASSERT_RELATION( BOX1::DIMENSION, ==, BOX2::DIMENSION );
    return MULTI_INDEX_BOX< BOX1::DIMENSION >(
        BOX1::MULTI_INDEX_TYPE::Componentwise_Max(box1.Min_Multi_Index(), box2.Min_Multi_Index()),
        BOX2::MULTI_INDEX_TYPE::Componentwise_Min(box1.Max_Multi_Index(), box2.Max_Multi_Index())
    );
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOX_INTERSECT_HPP
