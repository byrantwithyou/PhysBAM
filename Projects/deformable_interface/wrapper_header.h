#include "fixed_matrix.h"
#include "fixed_vector.h"
#include "polygon_mesh.h"

namespace data_exchange{
struct simulation_object
{
    int id;

    simulation_object(int i): id(i) {}
};

struct deformable_body : public simulation_object
{
    polygon_mesh mesh;

    std::vector<vf3> position;
    std::vector<vf3> velocity;
    float mass;

    static int fixed_id(int s = -1){static int i = s; return i;}
    deformable_body(): simulation_object(fixed_id()), mass(1) {}
};

struct ground_plane : public simulation_object
{
    vf3 position;
    vf3 normal;

    static int fixed_id(int s = -1){static int i = s; return i;}
    ground_plane(): simulation_object(fixed_id()), position(0,0,0), normal(0,1,0) {}
};

struct scripted_geometry : public simulation_object
{
    // Reference geometry
    polygon_mesh mesh;
    std::vector<vf3> position;

    static int fixed_id(int s = -1){static int i = s; return i;}
    scripted_geometry(): simulation_object(fixed_id()) {}
};

struct force
{
    int id;

    force(int i): id(i) {}
};

struct gravity_force : public force
{
    float magnitude;
    vf3 direction;

    static int fixed_id(int s = -1){static int i = s; return i;}
    gravity_force(): force(fixed_id()), magnitude(9.8f), direction(0,-1,0) {}
};

struct volumetric_force : public force
{
    float stiffness;
    float poissons_ratio;
    float damping;

    static int fixed_id(int s = -1){static int i = s; return i;}
    volumetric_force(): force(fixed_id()), stiffness(1e3f), poissons_ratio(.45f), damping(1.f) {}
};

inline void register_ids()
{
    int next_id = 0;
    deformable_body::fixed_id(next_id++);
    ground_plane::fixed_id(next_id++);
    scripted_geometry::fixed_id(next_id++);
    gravity_force::fixed_id(next_id++);
    volumetric_force::fixed_id(next_id++);
}
}


