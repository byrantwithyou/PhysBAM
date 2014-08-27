//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COUNT_BITS
//#####################################################################
#ifndef __COUNT_BITS__
#define __COUNT_BITS__

namespace PhysBAM{
inline int count_bits(unsigned d) {return d==0?0:1+count_bits(d&(d-1));}
}
#endif
