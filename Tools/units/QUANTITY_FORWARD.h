//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header QUANTITY_FORWARD
//#####################################################################
#ifndef __QUANTITY_FORWARD__
#define __QUANTITY_FORWARD__

#include <boost/mpl/bool.hpp>
#include <Tools/units/QUANTITY_FORWARD.h>
namespace PhysBAM{

namespace UNITS{
template<class T> class QUANTITY;
template<class T> struct CASTER;
}

using UNITS::QUANTITY;
using UNITS::CASTER;
namespace mpl=boost::mpl;

template<class T> struct IS_FLOAT_OR_DOUBLE;
template<class T> struct HAS_CHEAP_COPY;

template<class T> struct IS_FLOAT_OR_DOUBLE<QUANTITY<T> >:public mpl::true_{};
template<class T> struct HAS_CHEAP_COPY<QUANTITY<T> >:public mpl::true_{};

}
#endif
