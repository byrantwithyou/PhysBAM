//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "DEFORMABLE_EXAMPLE.h"
#include "simulation_object.h"
using namespace PhysBAM;
using namespace data_exchange;
#include<stdio.h>

struct physbam_simulation;
struct physbam_object;
struct physbam_force;

physbam_simulation * create_simulation()
{
    DEFORMABLE_EXAMPLE<float> * de = new DEFORMABLE_EXAMPLE<float>(STREAM_TYPE(0.f));
    de->Initialize_Simulation();
    return (physbam_simulation*) de;
}

void destroy_simulation(physbam_simulation * sim)
{
    DEFORMABLE_EXAMPLE<float> * de = (DEFORMABLE_EXAMPLE<float>*) sim;
    delete de;
}

physbam_object * add_object(physbam_simulation * sim, const simulation_object* so)
{
    DEFORMABLE_EXAMPLE<float> * de = (DEFORMABLE_EXAMPLE<float>*) sim;
    OBJECT_WRAPPER* wrap=de->Add_Simulation_Object(*so);
    return (physbam_object*) wrap;
}

physbam_force * add_force(physbam_simulation * sim, const force* f)
{
    DEFORMABLE_EXAMPLE<float> * de = (DEFORMABLE_EXAMPLE<float>*) sim;
    FORCE_WRAPPER* wrap=de->Add_Force(f);
    return (physbam_force*) wrap;
}

void apply_force_to_object(physbam_object * obj, physbam_force* f)
{
    OBJECT_WRAPPER * ow = (OBJECT_WRAPPER *) obj;
    FORCE_WRAPPER * fw = (FORCE_WRAPPER *) f;
    ow->de.new_forces_relations.Append(PAIR<OBJECT_WRAPPER*,FORCE_WRAPPER*>(ow,fw));
}

void simulate_frame(physbam_object * sim)
{
    DEFORMABLE_EXAMPLE<float> * de = (DEFORMABLE_EXAMPLE<float>*) sim;
    de->Simulate_Frame();
}



