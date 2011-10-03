#include <fstream>
#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "deformable_body_simulation.h"

// Compilation command: g++ -lboost_serialization main.cpp -g -o sample_test

using namespace data_exchange;
int main() {

    // Allocate the object that holds all of the data requred for the
    // simulation.
    deformable_body_simulation dbs;

    // To add a deformable body, allocate a deformable_body, populate it, and
    // add it to the simulation_objects array.  deformable_body_simulation will
    // free up the memory, so make sure you allocate the deformable_body with
    // new.  Polygons are specified in counterclockwise order.  This deformable
    // body has eight vertices and six faces.  It is a cube with square faces.
    deformable_body* db = new deformable_body;
    db->position.push_back(vf3(0,0,0));
    db->position.push_back(vf3(0,0,1));
    db->position.push_back(vf3(0,1,0));
    db->position.push_back(vf3(0,1,1));
    db->position.push_back(vf3(1,0,0));
    db->position.push_back(vf3(1,0,1));
    db->position.push_back(vf3(1,1,0));
    db->position.push_back(vf3(1,1,1));
    db->mesh.insert_polygon(vi4(7,6,2,3));
    db->mesh.insert_polygon(vi4(2,6,4,0));
    db->mesh.insert_polygon(vi4(1,0,4,5));
    db->mesh.insert_polygon(vi4(0,1,3,2));
    db->mesh.insert_polygon(vi4(3,1,5,7));
    db->mesh.insert_polygon(vi4(4,6,7,5));
    dbs.simulation_objects.push_back(db);

    // The next deformable body has six vertices and eight triangular faces.  It
    // is an octahedron.
    db = new deformable_body;
    db->position.clear();
    db->position.push_back(vf3(1,4,0));
    db->position.push_back(vf3(-1,4,0));
    db->position.push_back(vf3(0,5,0));
    db->position.push_back(vf3(0,3,0));
    db->position.push_back(vf3(0,4,1));
    db->position.push_back(vf3(0,4,-1));
    db->mesh.insert_polygon(vi3(0,2,4));
    db->mesh.insert_polygon(vi3(4,2,1));
    db->mesh.insert_polygon(vi3(0,5,2));
    db->mesh.insert_polygon(vi3(4,3,0));
    db->mesh.insert_polygon(vi3(1,3,4));
    db->mesh.insert_polygon(vi3(5,0,3));
    db->mesh.insert_polygon(vi3(2,5,1));
    db->mesh.insert_polygon(vi3(1,5,3));
    dbs.simulation_objects.push_back(db);

    // Next, add a ground.  The ground is slightly below the origin, and the
    // "up" direction is chosen to be 0,1,0, or the positive y axis.
    ground_plane* gp = new ground_plane;
    gp->position = vf3(0,-1,0);
    gp->normal = vf3(0,1,0);
    dbs.simulation_objects.push_back(gp);

    // This is almost like a deformable_object, except that it is only a surface
    // mesh.  In this case, we just make a copy of the geometry from the last
    // deformable body, but offset 3 units in the positive x direction.
    scripted_geometry* sc = new scripted_geometry;
    sc->mesh=db->mesh;
    sc->position=db->position;
    for(int i=0; i<sc->position.size(); i++) sc->position[i].data[0] += 3;
    dbs.simulation_objects.push_back(sc);

    // Now, we add a couple forces.  This force applies to the first two
    // objects, which are the two deformable bodies.  This gives the bodies
    // strength.
    volumetric_force* vf = new volumetric_force;
    vf->bodies_affected.push_back(0);
    vf->bodies_affected.push_back(1);
    vf->stiffness=1e2;
    dbs.simulation_forces.push_back(vf);

    // This forces makes the bodies fall.  It is applied to both deformable
    // bodies.
    gravity_force* gf = new gravity_force;
    gf->bodies_affected.push_back(0);
    gf->bodies_affected.push_back(1);
    dbs.simulation_forces.push_back(gf);

    // Write out the simulation object
    {
        std::ofstream ofs("filename");
        boost::archive::binary_oarchive oa(ofs);
        oa << dbs;
    }

    // Read it back in
    deformable_body_simulation dbs2;

    {
        std::ifstream ifs("filename");
        boost::archive::binary_iarchive ia(ifs);
        ia >> dbs2;
    }

    // Now, write it out a second time.  The two output files should match identically.
    {
        std::ofstream ofs2("filename2");
        boost::archive::binary_oarchive oa2(ofs2);
        oa2 << dbs2;
    }

    return 0;
}

