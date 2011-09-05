//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Unnur Gretarsdottir.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SKIN_DEFORMABLE_OBJECT_3D
//#####################################################################
#ifndef __SKIN_DEFORMABLE_OBJECT_3D__
#define __SKIN_DEFORMABLE_OBJECT_3D__

#include <Deformable_Objects/DEFORMABLE_OBJECT_3D.h>

namespace PhysBAM{
template<class T>
class SKIN_DEFORMABLE_OBJECT_3D:public DEFORMABLE_OBJECT_3D<T> {
public:
bool do_regular_collisions;
bool active;
SOLIDS_FLUIDS_EXAMPLE_3D<float, float> *example;

SKIN_DEFORMABLE_OBJECT_3D(): do_regular_collisions(true), active(true) 
    {}

    ~SKIN_DEFORMABLE_OBJECT_3D()
    {}


//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
void Advance_One_Time_Step(const T time,const T dt,const T cg_tolerance,const int cg_iterations,const bool perform_collision_body_collisions,const bool verbose,const int elastic_substeps=1)
{   
    std::cout<<"OK!!!!!!!!!!!!!!!!!!"<<std::endl;
    if(!active) return;
    // save velocity for later trapezoidal rule
    collisions.Reset_Object_Collisions();  //moved to beginning so data structures are available to postprocess_substep
    Save_Velocity(); 
    // update V implicitly to time n+1/2
    if(!Backward_Euler_Step_Velocity(dt/2,time,cg_tolerance,cg_iterations)){
        if(verbose) std::cout << "BACKWARD EULER FAILED!\n";Predictor_Corrector_Integrate_Velocity(time,time+dt/2);}
    // update position
    T dt_elastic=dt/elastic_substeps;
    Euler_Step_Position(dt_elastic);
    for(int i=2;i<=elastic_substeps;i++){Update_Position_Based_State();Euler_Step_Position(dt_elastic);}
    Restore_Velocity();
    // collision body collisions
    if(do_regular_collisions) { 
        if(perform_collision_body_collisions){
            external_forces_and_velocities->Activate_Collisions();
            int interactions=Adjust_Nodes_For_Collision_Body_Collisions(dt); 
            if(verbose && interactions) std::cout << "collision body collisions = " << interactions << std::endl;}
    } else {
        // In the inner layer sim, you can Snap to the Levelset either at the beginning
        // of each frame, or at each timestep, or both. If you want to do it at each
        // timestep, you should do it here...
        ((INNER_LAYER_SIM_EXAMPLE<T, float> *)example)->Snap_To_Levelset();
    }
    // positions just changed, so update again
    external_forces_and_velocities->Update_Time_Varying_Material_Properties(time,id_number);
    Update_Position_Based_State();
    // finish velocity update
    Euler_Step_Velocity(dt/2,time);
    if(!Backward_Euler_Step_Velocity(dt/2,time+dt/2,cg_tolerance,cg_iterations)){
        if(verbose) std::cout << "BACKWARD EULER FAILED!\n";Predictor_Corrector_Integrate_Velocity(time+dt/2,time+dt);}
    // reset collision velocities to zero
    external_forces_and_velocities->Activate_Collisions(false);
    // output
    if(verbose) std::cout << "maximum velocity  = " << particles.Maximum_Speed() << std::endl;
}
};
}
#endif
