//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_TAGGED_INDEX_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_TAGGED_INDEX_HPP

namespace PhysBAM
{

template< class T_TAG, class T_INDEX >
struct TAGGED_INDEX
{
    typedef T_TAG TAG_TYPE;
    typedef T_INDEX INDEX_TYPE;

    INDEX_TYPE value;

    static TAGGED_INDEX _(const INDEX_TYPE& value)
    {
        const TAGGED_INDEX tagged_index = { value };
        return tagged_index;
    }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_TAGGED_INDEX_HPP
