//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_REMOVE_IF_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_REMOVE_IF_HPP

#include <PhysBAM_Tools/Arrays/ARRAY.h>

namespace PhysBAM
{

template< class T, class T_PRED >
void
Remove_If(ARRAY<T>& a, T_PRED pred)
{
    for(int i = 1; i <= a.Size(); ++i) {
        if(!pred(a(i)))
            continue;
        int new_size = i - 1;
        while(++i <= a.Size())
            if(!pred(a(i)))
                a(++new_size) = a(i);
        a.Resize(new_size);
        return;
    }
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_REMOVE_IF_HPP
