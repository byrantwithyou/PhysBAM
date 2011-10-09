#ifndef __libmain__
#define __libmain__
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

extern "C"{
struct physbam_base {};
struct physbam_simulation : public physbam_base {};
struct physbam_object : public physbam_base {};
struct physbam_force : public physbam_base {};

physbam_simulation * create_simulation();
bool destroy_simulation(physbam_simulation * sim);
physbam_object * add_object(physbam_simulation * sim, const data_exchange::simulation_object* so);
physbam_force * add_force(physbam_simulation * sim, const data_exchange::force* f);
bool apply_force_to_object(physbam_object * obj, physbam_force* f);
bool simulate_frame(physbam_simulation * sim);

void set_int(physbam_base * obj, const char * attribute, int x);
int get_int(const physbam_base * obj, const char * attribute);
void set_float(physbam_base * obj, const char * attribute, float x);
float get_float(const physbam_base * obj, const char * attribute);
int get_array_length(const physbam_base * obj, const char * attribute);
void set_int_array(physbam_base * obj, const char * attribute, const int * x, int length, int start = 0);
void get_int_array(const physbam_base * obj, const char * attribute, int * x, int length, int start = 0);
void set_float_array(physbam_base * obj, const char * attribute, const float * x, int length, int start = 0);
void get_float_array(const physbam_base * obj, const char * attribute, float * x, int length, int start = 0);
}

#endif
