//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KINEMATIC_EVOLUTION
//#####################################################################
#ifndef __KINEMATIC_EVOLUTION__
#define __KINEMATIC_EVOLUTION__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Matrices/MATRIX_POLICY.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGID_BODY_STATE;
template<class TV> class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES;

template<class TV>
class KINEMATIC_EVOLUTION
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
public:
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>& rigid_body_example_velocities;

    KINEMATIC_EVOLUTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>& rigid_body_example_velocities);
    virtual ~KINEMATIC_EVOLUTION();

//#####################################################################
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time);
    virtual void Set_External_Velocities(TWIST<TV>& twist,const T time,const int id);
    virtual void Set_Kinematic_Velocities(TWIST<TV>& twist,const T frame_dt,const T time,const int id);
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time);
    virtual void Set_External_Positions(FRAME<TV>& frame,const T time,const int id);
    void Reset_Kinematic_Rigid_Bodies(const T time);
//#####################################################################
};
}
#endif
