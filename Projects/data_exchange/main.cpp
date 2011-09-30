#include <fstream>
#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "deformable_body_simulation.h"

using namespace data_exchange;
int main() {
    deformable_body_simulation dbs;

    deformable_body* db = new deformable_body;
    ground_plane* gp = new ground_plane;
    scripted_geometry* sc = new scripted_geometry;

    dbs.simulation_objects.push_back(db);
    dbs.simulation_objects.push_back(gp);
    dbs.simulation_objects.push_back(sc);

    volumetric_force* vf = new volumetric_force;
    gravity_force* gf = new gravity_force;

    dbs.simulation_forces.push_back(vf);
    dbs.simulation_forces.push_back(gf);

    // create and open a character archive for output
    std::ofstream ofs("filename");
    std::ifstream ifs("filename");
    std::ofstream ofs2("filename2");
    boost::archive::text_oarchive oa(ofs);
    boost::archive::text_iarchive ia(ifs);
    boost::archive::text_oarchive oa2(ofs2);

    oa << dbs;

    deformable_body_simulation dbs2;

    ia >> dbs2;

    oa2 << dbs2;

    return 0;
}

