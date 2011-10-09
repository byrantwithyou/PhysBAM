#include "deformable_body.h"
#include "gravity_force.h"
#include "ground_plane.h"
#include "scripted_geometry.h"
#include "simulation_object.h"
#include "volumetric_force.h"
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>

struct physbam_simulation;
struct physbam_object;
struct physbam_force;

int main(int argc,char* argv[])
{
    using namespace data_exchange;
    void* handle = dlopen("libPhysBAM_Wrapper.so", RTLD_LAZY);
    if(!handle){
        const char *p = dlerror();
        printf("dlopen error %s\n", p);
        exit(1);
    }

    physbam_simulation * (*create_simulation)() = (physbam_simulation* (*)()) dlsym(handle, "create_simulation");
    bool (*destroy_simulation)(physbam_simulation*) = (bool (*)(physbam_simulation*)) dlsym(handle, "destroy_simulation");
    physbam_object* (*add_object)(physbam_simulation*, const simulation_object*) = (physbam_object* (*)(physbam_simulation*, const simulation_object*)) dlsym(handle, "add_object");
    physbam_force* (*add_force)(physbam_simulation*, const force*) = (physbam_force* (*)(physbam_simulation*, const force*)) dlsym(handle, "add_force");
    bool (*apply_force_to_object)(physbam_object*, physbam_force*) = (bool (*)(physbam_object*, physbam_force*)) dlsym(handle, "apply_force_to_object");
    bool (*simulate_frame)(physbam_simulation*) = (bool (*)(physbam_simulation*)) dlsym(handle, "simulate_frame");

    printf("%p %p %p %p %p %p\n", create_simulation, destroy_simulation, add_object, add_force, apply_force_to_object, simulate_frame);

    physbam_simulation * sim = create_simulation();

    deformable_body db;
    db.position.push_back(vf3(0,0,0));
    db.position.push_back(vf3(0,0,1));
    db.position.push_back(vf3(0,1,0));
    db.position.push_back(vf3(0,1,1));
    db.position.push_back(vf3(1,0,0));
    db.position.push_back(vf3(1,0,1));
    db.position.push_back(vf3(1,1,0));
    db.position.push_back(vf3(1,1,1));
    db.mesh.insert_polygon(vi4(7,6,2,3));
    db.mesh.insert_polygon(vi4(2,6,4,0));
    db.mesh.insert_polygon(vi4(1,0,4,5));
    db.mesh.insert_polygon(vi4(0,1,3,2));
    db.mesh.insert_polygon(vi4(3,1,5,7));
    db.mesh.insert_polygon(vi4(4,6,7,5));
    physbam_object* d1 = add_object(sim, &db);

    db.position.clear();
    db.position.push_back(vf3(1,4,0));
    db.position.push_back(vf3(-1,4,0));
    db.position.push_back(vf3(0,5,0));
    db.position.push_back(vf3(0,3,0));
    db.position.push_back(vf3(0,4,1));
    db.position.push_back(vf3(0,4,-1));
    db.mesh.insert_polygon(vi3(0,2,4));
    db.mesh.insert_polygon(vi3(4,2,1));
    db.mesh.insert_polygon(vi3(0,5,2));
    db.mesh.insert_polygon(vi3(4,3,0));
    db.mesh.insert_polygon(vi3(1,3,4));
    db.mesh.insert_polygon(vi3(5,0,3));
    db.mesh.insert_polygon(vi3(2,5,1));
    db.mesh.insert_polygon(vi3(1,5,3));
    physbam_object* d2 = add_object(sim, &db);

    ground_plane gp;
    gp.position = vf3(0,-1,0);
    gp.normal = vf3(0,1,0);
    add_object(sim, &gp);

    scripted_geometry sc;
    sc.mesh=db.mesh;
    sc.position=db.position;
    for(size_t i=0; i<sc.position.size(); i++) sc.position[i].data[0] += 3;
    add_object(sim, &sc);

    volumetric_force vf;
    vf.stiffness=1e3;
    vf.damping=.01;
    physbam_force* f1 = add_force(sim, &vf);

    gravity_force gf;
    physbam_force* f2 = add_force(sim, &gf);

    apply_force_to_object(d1, f1);
    apply_force_to_object(d2, f1);
    apply_force_to_object(d1, f2);
    apply_force_to_object(d2, f2);

    simulate_frame(sim);
    simulate_frame(sim);
    simulate_frame(sim);
    simulate_frame(sim);
    simulate_frame(sim);

    destroy_simulation(sim);

    return 0;
}
