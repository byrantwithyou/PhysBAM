#include <fstream>
#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "deformable_body_simulation.h"

// (setq compile-command "g++ -lboost_serialization main.cpp -g")

using namespace data_exchange;
int main() {
    deformable_body_simulation dbs;

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

    ground_plane* gp = new ground_plane;
    gp->position = vf3(0,-1,0);
    gp->normal = vf3(0,1,0);
    dbs.simulation_objects.push_back(gp);

    scripted_geometry* sc = new scripted_geometry;
    sc->mesh=db->mesh;
    sc->position=db->position;
    for(int i=0; i<sc->position.size(); i++) sc->position[i].data[0] += 2;
    dbs.simulation_objects.push_back(sc);

    volumetric_force* vf = new volumetric_force;
    dbs.simulation_forces.push_back(vf);

    gravity_force* gf = new gravity_force;
    dbs.simulation_forces.push_back(gf);

    {
        std::ofstream ofs("filename");
        boost::archive::text_oarchive oa(ofs);
        oa << dbs;
    }

    deformable_body_simulation dbs2;

    {
        std::ifstream ifs("filename");
        boost::archive::text_iarchive ia(ifs);
        ia >> dbs2;
    }

    {
        std::ofstream ofs2("filename2");
        boost::archive::text_oarchive oa2(ofs2);
        oa2 << dbs2;
    }

    return 0;
}

