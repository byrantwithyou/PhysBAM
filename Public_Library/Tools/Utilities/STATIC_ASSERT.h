//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header STATIC_ASSERT
//#####################################################################
#ifndef __STATIC_ASSERT__
#define __STATIC_ASSERT__

namespace PhysBAM{

#define STATIC_ASSERT(...) static_assert((__VA_ARGS__),#__VA_ARGS__)
}

#endif
