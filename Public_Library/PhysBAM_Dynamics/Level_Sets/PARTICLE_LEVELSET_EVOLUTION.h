//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_EVOLUTION
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_EVOLUTION__
#define __PARTICLE_LEVELSET_EVOLUTION__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <cassert>
namespace PhysBAM{

template<class T>
class PARTICLE_LEVELSET_EVOLUTION:public NONCOPYABLE
{
public:
    bool use_particle_levelset;
    T time;
    T cfl_number;
    bool use_reinitialization,use_fmm;
    int reseeding_frequency;
    bool use_frozen_velocity;
    int runge_kutta_order_particles,runge_kutta_order_levelset;
    bool track_mass;
    double initial_mass;

    PARTICLE_LEVELSET_EVOLUTION()
        :use_particle_levelset(true),time(0),cfl_number((T).5),use_reinitialization(false),use_fmm(true),
        reseeding_frequency(20),use_frozen_velocity(true),runge_kutta_order_particles(2),
        runge_kutta_order_levelset(1),track_mass(),initial_mass()
    {}

    virtual ~PARTICLE_LEVELSET_EVOLUTION()
    {}

    virtual void Set_CFL_Number(const T cfl_number_input)
    {cfl_number=cfl_number_input;}

    void Use_Reinitialization()
    {use_reinitialization=true;use_fmm=false;}

    void Use_Fast_Marching_Method()
    {use_fmm=true;use_reinitialization=false;}

    void Set_Runge_Kutta_Order_Levelset(const int order)
    {runge_kutta_order_levelset=order;assert(order>=1 && order<=3);}
   
    void Set_Runge_Kutta_Order_Particles(const int order)
    {runge_kutta_order_particles=order;assert(order>=1 && order<=3);}

//#####################################################################
};
}
#endif
