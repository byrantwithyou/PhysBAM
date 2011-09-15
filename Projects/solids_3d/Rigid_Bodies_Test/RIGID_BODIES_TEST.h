//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODIES_TEST
//#####################################################################
#ifndef __RIGID_BODIES_TEST__
#define __RIGID_BODIES_TEST__

#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template <class T,class RW>
class RIGID_BODIES_TEST:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart_frame;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::frame_rate;

    RIGID_BODIES_TEST()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,fluids_parameters.NONE)
    {
        last_frame=10*24;
        restart=false;restart_frame=0;   
        frame_rate=24;
        solids_parameters.cfl=(T).5;
        output_directory="Rigid_Bodies_Test/output";
        verbose_dt=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.gravity=9.8;
        solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=1;
    }

    ~RIGID_BODIES_TEST()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int id;
#if 0
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere");
//    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(0,10,0);
//    solids_parameters.rigid_body_parameters.list(id)->twist.linear=VECTOR<T,3>(2,0,0);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(0);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(0);
#endif

#if 0
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(5,10,0);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(0.5);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(0);
#endif

#if 0
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(10,10,0);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(1);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(0);
#endif

#if 1
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/subdivided_box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(0,6,0);
////    solids_parameters.rigid_body_parameters.list(id)->frame.r=QUATERNION<T>(0.01,VECTOR<T,3>(1,1,1));
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(0);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(1);
#endif

#if 1
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/subdivided_box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
    solids_parameters.rigid_body_parameters.list(id)->twist.linear=VECTOR<T,3>(0,6,0);
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(0,3,0);
    solids_parameters.rigid_body_parameters.list(id)->frame.r=QUATERNION<T>(0.1,VECTOR<T,3>(0,1,0));
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(0);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(1);
#endif

#if 0
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/subdivided_box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(25,10,0);
//    solids_parameters.rigid_body_parameters.list(id)->frame.r=QUATERNION<T>(0.01,VECTOR<T,3>(1,1,1));
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(1);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(1);
#endif

#if 0
    for(int i=1;i<10;i++){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/subdivided_box");
        solids_parameters.rigid_body_parameters.list(id)->Set_Name("second object");
        solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(3*i,10+0.1*i,0+0.1*i);
        solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    }
#endif

#if 0
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("second object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(5,10.1,0.1);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("third object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(10,10.2,0.2);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("fourth object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,3>(15,10.3,0.3);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
#endif

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T)1);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(1);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
#if 1
    solids_evolution->rigid_body_collisions->collision_manager.Use_Collision_Matrix();
//    solids_evolution->rigid_body_collisions->collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(1,2,false);
    solids_evolution->rigid_body_collisions->collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(2,1,false);
//    solids_evolution->rigid_body_collisions->collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(3,4,false);
#endif
}
//#####################################################################
};
}
#endif
