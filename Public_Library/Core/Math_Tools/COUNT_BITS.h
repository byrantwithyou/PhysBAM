//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COUNT_BITS
//#####################################################################
#ifndef __COUNT_BITS__
#define __COUNT_BITS__

namespace PhysBAM{
constexpr inline int count_bits(unsigned a)
{
    a = (a & 0x55555555) + ((a>>1) & 0x55555555);
    a = (a & 0x33333333) + ((a>>2) & 0x33333333);
    a = (a & 0x0f0f0f0f) + ((a>>4) & 0x0f0f0f0f);
    a = (a & 0x00ff00ff) + ((a>>8) & 0x00ff00ff);
    a = (a & 0x0000ffff) + ((a>>16) & 0x0000ffff);
    return a;
}
constexpr inline int count_bits(int a){return count_bits((unsigned)a);}
constexpr inline int count_bits(unsigned long long a)
{
    a = (a & 0x5555555555555555) + ((a>>1) & 0x5555555555555555);
    a = (a & 0x3333333333333333) + ((a>>2) & 0x3333333333333333);
    a = (a & 0x0f0f0f0f0f0f0f0f) + ((a>>4) & 0x0f0f0f0f0f0f0f0f);
    a = (a & 0x00ff00ff00ff00ff) + ((a>>8) & 0x00ff00ff00ff00ff);
    a = (a & 0x0000ffff0000ffff) + ((a>>16) & 0x0000ffff0000ffff);
    a = (a & 0x00000000ffffffff) + ((a>>32) & 0x00000000ffffffff);
    return a;
}
}
#endif
