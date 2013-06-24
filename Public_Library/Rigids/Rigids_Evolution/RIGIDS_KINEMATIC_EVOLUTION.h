//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_KINEMATIC_EVOLUTION
//#####################################################################
#ifndef __RIGIDS_KINEMATIC_EVOLUTION__
#define __RIGIDS_KINEMATIC_EVOLUTION__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Matrices/MATRIX_POLICY.h>
#include <Geometry/Solids_Geometry_Evolution/KINEMATIC_EVOLUTION.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class RIGID_BODY_STATE;

template<class TV>
class RIGIDS_KINEMATIC_EVOLUTION:public KINEMATIC_EVOLUTION<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
    typedef KINEMATIC_EVOLUTION<TV> BASE;
public:
    using BASE::Set_External_Velocities;using BASE::Set_Kinematic_Velocities;using BASE::Set_External_Positions;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    RIGIDS_KINEMATIC_EVOLUTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>& rigid_body_example_velocities,bool use_kinematic_keyframes_input);
    ~RIGIDS_KINEMATIC_EVOLUTION();

//#####################################################################
    void Set_External_Velocities(TWIST<TV>& twist,const T time,const int id);
    void Set_Kinematic_Velocities(TWIST<TV>& twist,const T frame_dt,const T time,const int id);
    void Set_External_Positions(FRAME<TV>& frame,const T time,const int id);
//#####################################################################
};
}
#endif
