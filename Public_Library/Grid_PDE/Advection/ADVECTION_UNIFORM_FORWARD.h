//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_UNIFORM_FORWARD
//#####################################################################
#ifndef __ADVECTION_UNIFORM_FORWARD__
#define __ADVECTION_UNIFORM_FORWARD__

#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM_FORWARD.h>
namespace PhysBAM{

template<class TV,class T2,class T_AVERAGING=AVERAGING_UNIFORM<TV>,class T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<TV,T2> > class ADVECTION_SEMI_LAGRANGIAN_UNIFORM;
template<class TV,class T2,class T_AVERAGING=AVERAGING_UNIFORM<TV>,class T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<TV,T2> > class ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA;

template<class TV,class T2,class T_NESTED_ADVECTION> class ADVECTION_MACCORMACK_UNIFORM;

template<class TV,class T2,class T_AVERAGING=AVERAGING_UNIFORM<TV> > class ADVECTION_SEPARABLE_UNIFORM;
template<class T,class T2> class ADVECTION_CONSERVATIVE_ENO;
template<class T,class T2> class ADVECTION_CONSERVATIVE_WENO;
template<class TV,class T2> class ADVECTION_CENTRAL;
template<class TV,class T2> class ADVECTION_HAMILTON_JACOBI_ENO;
template<class TV,class T2> class ADVECTION_HAMILTON_JACOBI_WENO;
template<class T,class T2> class ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO;
template<class T,class T2> class ADVECTION_NONCONSERVATIVE_FLUX_FORM_WENO;

}
#endif
