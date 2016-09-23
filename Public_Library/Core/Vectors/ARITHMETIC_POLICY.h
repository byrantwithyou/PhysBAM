//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARITHMETIC_POLICY
//#####################################################################
#ifndef __ARITHMETIC_POLICY__
#define __ARITHMETIC_POLICY__

#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Core/Vectors/SCALAR_POLICY.h>
#include <Core/Vectors/VECTOR_FORWARD.h>

#ifdef DIFFERENCE
#undef DIFFERENCE
#endif

namespace PhysBAM{

template<class T1,class T2,class ENABLER=void> struct CAN_ASSIGN;

// Builtin
template<class T> struct CAN_ASSIGN<T,T>{static const bool value=true;};
template<> struct CAN_ASSIGN<double,float>{static const bool value=false;};
template<> struct CAN_ASSIGN<double,int>{static const bool value=false;};
template<> struct CAN_ASSIGN<float,double>{static const bool value=false;};
template<> struct CAN_ASSIGN<float,int>{static const bool value=false;};
template<> struct CAN_ASSIGN<int,double>{static const bool value=false;};
template<> struct CAN_ASSIGN<int,float>{static const bool value=false;};

}
#endif
