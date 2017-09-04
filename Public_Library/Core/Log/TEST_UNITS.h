//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TEST_UNITS__
#define __TEST_UNITS__
#include <Core/Log/LOG.h>
namespace PhysBAM{

#define TEST_UNITS(x) LOG::printf("TEST_UNITS %s:%i %s %s @BEGIN@ %.16P @END@\n",__FILE__,__LINE__,__FUNCTION__,#x,x)

}
#endif
