//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LAPLACE_POLICY__
#define __LAPLACE_POLICY__

namespace PhysBAM{
template<class TV> class LAPLACE_UNIFORM;

template<class TV>
struct LAPLACE_POLICY
{
    typedef LAPLACE_UNIFORM<TV> LAPLACE;
};
}
#endif
