//#####################################################################
// Copyright 2004, Eran_Guendelman, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KINEMATIC_EXAMPLE
//##################################################################### 
#ifndef __KINEMATIC_EXAMPLE__
#define __KINEMATIC_EXAMPLE__

namespace PhysBAM{

template<class T,class RW>
class KINEMATIC_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::output_directory;using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;

    T time1,time2,time3,time4,time5;
    T kinematic_body_vel;

    KINEMATIC_EXAMPLE()
        :BASE(0,fluids_parameters.NONE)
    {
        last_frame=600;
        output_directory="Kinematic_Example/output";
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    T baseboxsize=2;
    T dropboxsize=(T)1.2;
    T dropheight=70;
    T smallboxsize=(T)0.5;
    T smallboxmass=1;
    T offset=(T)0.05;

    T stack_epsilon=(T)0.3;
    T stack_mu=(T)0.5;

    time1=(T)10/24;
    time2=(T)200/24;
    time3=(T)350/24;
    time4=(T)400/24;

    T boxsize1=(T)1.0;
    T boxsize2=(T)0.5;

    int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",boxsize1);
    solids_parameters.rigid_body_parameters.list(index)->frame.t=VECTOR_3D<T>(0,2*baseboxsize+boxsize1,0);
    solids_parameters.rigid_body_parameters.list(index)->Set_Coefficient_Of_Restitution(stack_epsilon);
    solids_parameters.rigid_body_parameters.list(index)->Set_Coefficient_Of_Friction(stack_mu);
    solids_parameters.rigid_body_parameters.list(index)->Set_Mass(smallboxmass);
    solids_parameters.rigid_body_parameters.list(index)->Set_Name("stack box 1a");
    solids_parameters.rigid_body_parameters.list(index)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",boxsize2);
    solids_parameters.rigid_body_parameters.list(index)->frame.t=VECTOR_3D<T>(0,2*baseboxsize+2*boxsize1+boxsize2,0);
    solids_parameters.rigid_body_parameters.list(index)->Set_Coefficient_Of_Restitution(stack_epsilon);
    solids_parameters.rigid_body_parameters.list(index)->Set_Coefficient_Of_Friction(stack_mu);
    solids_parameters.rigid_body_parameters.list(index)->Set_Mass(smallboxmass);
    solids_parameters.rigid_body_parameters.list(index)->Set_Name("stack box 1b");
    solids_parameters.rigid_body_parameters.list(index)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",baseboxsize);
    solids_parameters.rigid_body_parameters.list(index)->Set_Coefficient_Of_Restitution((T)0.1);
    solids_parameters.rigid_body_parameters.list(index)->Set_Coefficient_Of_Friction(stack_mu);
    solids_parameters.rigid_body_parameters.list(index)->Set_Name("base box");
    solids_parameters.rigid_body_parameters.list(index)->is_kinematic = true;
    kinematic_body_vel=5;

    for(int i=1;i<=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;i++){RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(i);
        if(rigid_body.is_kinematic) Set_Kinematic_Positions(rigid_body.frame,time,rigid_body.particle_index);}
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    T baseboxsize=2;
    if(time<time1){frame.t=VECTOR_3D<T>(0,baseboxsize,0);frame.r=QUATERNION<T>();}
    else if(time<time2){frame.t=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time-time1),0);frame.r=QUATERNION<T>();}
    else if(time<time3){frame.t=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1),0);frame.r=QUATERNION<T>();}
    else if(time<time4){frame.t=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1),0);frame.r=QUATERNION<T>(-(T)0.4*(time-time3),VECTOR_3D<T>(0,0,1));}
    else{frame.t=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1),0);frame.r=QUATERNION<T>(-(T)0.4*(time4-time3),VECTOR_3D<T>(0,0,1));}
}
//#####################################################################
};
}
#endif
