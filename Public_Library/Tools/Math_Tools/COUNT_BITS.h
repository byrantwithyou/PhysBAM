//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COUNT_BITS
//#####################################################################
#ifndef __COUNT_BITS__
#define __COUNT_BITS__

namespace PhysBAM{
template<unsigned d> struct COUNT_BITS;

template<> struct COUNT_BITS<0>{enum WORKAROUND {value=0};};

template<unsigned d> struct COUNT_BITS
{
    enum WORKAROUND {value=1+COUNT_BITS<(d&(d-1))>::value};
};
}
#endif
