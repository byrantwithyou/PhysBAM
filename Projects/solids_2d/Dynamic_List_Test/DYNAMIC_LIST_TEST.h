//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYNAMIC_LIST_TEST
//#####################################################################
#ifndef __DYNAMIC_LIST_TEST__
#define __DYNAMIC_LIST_TEST__

#include <Rigid_Bodies/RIGID_BODY_IMPULSE_ACCUMULATOR_2D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template <class T,class RW>
class DYNAMIC_LIST_TEST:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::verbose_dt;

    DYNAMIC_LIST_TEST()
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.NONE)
    {
        last_frame=10*24;
        restart=false;restart_frame=0;   
        solids_parameters.cfl=(T).5;
        output_directory="Dynamic_List_Test/output";
        verbose_dt=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

    ~DYNAMIC_LIST_TEST()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies_2D/ground");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square");
    solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
    solids_parameters.rigid_body_parameters.list(id)->Set_Mass(1000);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>(-3,5);
    solids_parameters.rigid_body_parameters.list(id)->velocity=VECTOR_2D<T>(10,0);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Update_Animated_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE{
#if 1
    static int state=1;
    if (time > 2 && state==1){state++;
        std::cout << "Adding second body" << std::endl;
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square");
        solids_parameters.rigid_body_parameters.list(id)->Set_Name("second object");
        solids_parameters.rigid_body_parameters.list(id)->Set_Mass(1000);
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>(0,10);
        solids_parameters.rigid_body_parameters.list(id)->velocity=VECTOR_2D<T>(0,0);
        solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    }
    else if (time > 4 && state==2){state++;
        std::cout << "Adding third body" << std::endl;
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/circle");
        solids_parameters.rigid_body_parameters.list(id)->Set_Name("third object");
        solids_parameters.rigid_body_parameters.list(id)->Set_Mass(1000);
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>(0,10);
        solids_parameters.rigid_body_parameters.list(id)->velocity=VECTOR_2D<T>(0,0);
        solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    }
    else if (time > 6 && state==3){state++;
        std::cout << "Removing second body" << std::endl;
        solids_parameters.rigid_body_parameters.list.Remove_Element(2);
    }
    else if (time > 8 && state==4){state++;
        std::cout << "Adding fourth body" << std::endl;
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square");
        solids_parameters.rigid_body_parameters.list(id)->Set_Name("fourth object");
        solids_parameters.rigid_body_parameters.list(id)->Set_Mass(1000);
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>(0,10);
        solids_parameters.rigid_body_parameters.list(id)->velocity=VECTOR_2D<T>(0,0);
        solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    }
#endif
}
//#####################################################################
};
}
#endif
