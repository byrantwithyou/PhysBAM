//#####################################################################
// Copyright 2004, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KINEMATIC_EXAMPLE
//##################################################################### 
#ifndef __KINEMATIC_EXAMPLE__
#define __KINEMATIC_EXAMPLE__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include "../RIGID_BODIES_2D_EXAMPLE.h"
#include <Forces_And_Torques/EXTERNAL_FORCES_AND_VELOCITIES.h>
namespace PhysBAM{

template<class T>
class KINEMATIC_EXAMPLE:public RIGID_BODIES_2D_EXAMPLE<T>,public EXTERNAL_FORCES_AND_VELOCITIES<T,VECTOR_2D<T>,T,T>
{
public:
    T time1,time2,time3,time4,time5;
    RIGID_BODY_STATE_2D<T> current_state,next_state;
    int kinematic_body_index;
    T kinematic_body_vel;

    KINEMATIC_EXAMPLE(int parameter=3)
        :RIGID_BODIES_2D_EXAMPLE<T>(parameter)
    {
        last_frame=240;
        output_directory="Kinematic_Example/output";
        if(parameter==3) artificial_maximum_speed=20;
    }
//#####################################################################
    void Initialize_Rigid_Bodies() PHYSBAM_OVERRIDE;
    void Update_Animated_Parameters(const T time) PHYSBAM_OVERRIDE;
    void Compute_State(RIGID_BODY_STATE_2D<T>& state,const T time);
    void Compute_Velocities(const RIGID_BODY<TV>& rigid_body,RIGID_BODY_STATE_2D<T>& state1,RIGID_BODY_STATE_2D<T>& state2);
    void Set_External_Velocities(VECTOR_2D<T>& V,T& omega,const T time,const int id_number=1);
    void Set_External_Position_And_Orientation(VECTOR_2D<T>& X,T& orientation,const T time,const int id_number=1);
//#####################################################################
};
}
#endif
