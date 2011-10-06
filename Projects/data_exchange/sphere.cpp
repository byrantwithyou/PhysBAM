#include <fstream>
#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "deformable_body_simulation.h"

// Compilation command: g++ -lboost_serialization sphere.cpp -g -o sample_test

using namespace data_exchange;
int main() {

    deformable_body_simulation dbs;

    int n=40,m=20;
    deformable_body* db = new deformable_body;
    float pi=4*atan(1);

    for(int i=0; i<n; i++)
      {
        int ni=(i+1)%n;
        for(int j=0; j<m; j++)
          {
            int nj=(j+1)%m;
            db->position.push_back(vf3(cos(2*pi*i/n)*sin(pi*(j+1)/(m+1)),sin(2*pi*i/n)*sin(pi*(j+1)/(m+1)),cos(pi*(j+1)/(m+1))));
            if(j<m-1) db->mesh.insert_polygon(vi4(i*m+j,i*m+nj,ni*m+nj,ni*m+j));
          }
        db->mesh.insert_polygon(vi3(ni*m+m-1,i*m+m-1,m*n));
        db->mesh.insert_polygon(vi3(i*m,ni*m,m*n+1));
      }
    db->position.push_back(vf3(0,0,-1));
    db->position.push_back(vf3(0,0,1));
    for(int i=0; i<db->position.size(); i++) db->position[i].data[1] += 3;
    dbs.simulation_objects.push_back(db);

    ground_plane* gp = new ground_plane;
    gp->position = vf3(0,-1,0);
    gp->normal = vf3(0,1,0);
    dbs.simulation_objects.push_back(gp);

    scripted_geometry* sc = new scripted_geometry;
    sc->mesh=db->mesh;
    sc->position=db->position;
    for(int i=0; i<sc->position.size(); i++) sc->position[i].data[1] -= 3;
    for(int i=0; i<sc->position.size(); i++) sc->position[i].data[0] += .1;
    dbs.simulation_objects.push_back(sc);

    sc = new scripted_geometry;
    sc->mesh=db->mesh;
    sc->position=db->position;
    for(int i=0; i<sc->position.size(); i++) sc->position[i].data[1] -= 3;
    for(int i=0; i<sc->position.size(); i++) sc->position[i].data[0] -= 4.5;
    dbs.simulation_objects.push_back(sc);

    volumetric_force* vf = new volumetric_force;
    vf->bodies_affected.push_back(0);
    vf->stiffness=1e3;
    vf->damping=.01;
    dbs.simulation_forces.push_back(vf);

    gravity_force* gf = new gravity_force;
    gf->bodies_affected.push_back(0);
    dbs.simulation_forces.push_back(gf);

    {
        std::ofstream ofs("filename");
        boost::archive::binary_oarchive oa(ofs);
        oa << dbs;
    }

    return 0;
}

