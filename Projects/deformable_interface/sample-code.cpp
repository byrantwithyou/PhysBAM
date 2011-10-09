--------------------- PhysLibrary part -------------------------

physLib.h

#ifndef __PHYSLIB_H__
#define __PHYSLIB_H__

class PhysLib
{
public:

PhysLib();

/* use virtual otherwise linker will try to perform static linkage */

virtual void SetParameter(char *parm , int val);
    virtual int GetParameter(char *parm);

    virtual bool RunStep(float time_step);

private:
 int x;
};

physLib.cpp

#include "physLib.h"
#include <iostream>

extern "C" PhysLib* create_simulation(){
 return new PhysLib;
}

extern "C" void destroy_simualtion( PhysLib* object ){
 delete object;
}

PhysLib::PhysLib(){
 // init simulation here
}

void PhysLib::SetParameter(char *parm, int val){
    parms[parm] = val;
}

int PhysLib::GetParameters(char *parm){
    return parms[parm];
}

----------------------- Houdini Solver part ---------------------------

physSolver.cpp

#include <dlfcn.h>
#include <iostream>
#include "physLib.h"

// some code here ------

void* handle = dlopen("physlib.so", RTLD_LAZY);

physLib* (*create)();
void (*destroy)(physLib*);

create_sim = (physLib* (*)())dlsym( handle, "create_simulation" );
destroy_sim = (void (*)(physLib*))dlsym( handle, "destroy_simulation" );

physLib *sim = (physLib*)create_simulation();
physLib->SetParameter( “simulation_speed” , 1.0f );

physObj *obj = sim->AddObject( TRI_MESH );
obj->SetPoints( some_points_array );
obj->SetPolygons( some_point_indexes_array );
obj->SetParameter( “mass”, 100.0f );

physForce *force = sim->AddForce( SIM_GRAVITY );
force->SetParameter( “direction”, vec3(0,1,0) );

obj->AffectByForce( force );

sim->RunStep( 0.2f );

const physPointsArray *pts = obj->ReadBackPoints();
// update Houdini’s geometry here
// or read new object transform back in case of rigid body like:
const physTransform *xform = obj->ReadBackTransform();
// update houdini object transform here

destroy( sim );

// some code here

