//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODIES_TEST
//#####################################################################
#ifndef __RIGID_BODIES_TEST__
#define __RIGID_BODIES_TEST__

#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template <class T,class RW>
class RIGID_BODIES_TEST:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;
    using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;using BASE::verbose_dt;

    RIGID_BODIES_TEST()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,fluids_parameters.NONE)
    {
        last_frame=10*24;
        restart=false;restart_frame=0;   
        solids_parameters.cfl=(T).5;
        output_directory="Rigid_Bodies_Test/output";
        verbose_dt=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.rigid_body_parameters.use_collision_matrix=true;
        solids_parameters.gravity=0;
    }

    ~RIGID_BODIES_TEST()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,2>(0,10);
    solids_parameters.rigid_body_parameters.list(id)->twist.linear=VECTOR<T,2>(2,0);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("second object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,2>(5,10.1);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("third object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,2>(10,10.2);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("fourth object");
    solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR<T,2>(15,10.3);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies_2D/ground");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
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
    solids_evolution->rigid_body_collisions->collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(1,2,false);
    solids_evolution->rigid_body_collisions->collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(2,1,false);
    solids_evolution->rigid_body_collisions->collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(3,4,false);
}
//#####################################################################
};
}
#endif
