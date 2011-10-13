#ifndef __libmain__
#define __libmain__
#include "wrapper_header.h"

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

int get_id(physbam_base * obj, const char * attribute);
void set_int(physbam_base * obj, int id, int x);
int get_int(const physbam_base * obj, int id);
void set_float(physbam_base * obj, int id, float x);
float get_float(const physbam_base * obj, int id);
void set_vf3(physbam_base * obj, int id, data_exchange::vf3 x);
data_exchange::vf3 get_vf3(const physbam_base * obj, int id);
int get_array_length(const physbam_base * obj, int id);
int get_int_array_length(const physbam_base * obj, int id);
void set_int_array(physbam_base * obj, int id, const int * x, int length, int start);
void get_int_array(const physbam_base * obj, int id, int * x, int length, int start);
int get_float_array_length(const physbam_base * obj, int id);
void set_float_array(physbam_base * obj, int id, const float * x, int length, int start);
void get_float_array(const physbam_base * obj, int id, float * x, int length, int start);
int get_vf3_array_length(const physbam_base * obj, int id);
void set_vf3_array(physbam_base * obj, int id, const data_exchange::vf3 * x, int length, int start);
void get_vf3_array(const physbam_base * obj, int id, data_exchange::vf3 * x, int length, int start);
}

enum funcid_type
{
    funcid_create_simulation, funcid_destroy_simulation, funcid_add_object, funcid_add_force,
    funcid_apply_force_to_object, funcid_simulate_frame, funcid_get_id, funcid_set_int,
    funcid_get_int, funcid_set_float, funcid_get_float, funcid_set_vf3,
    funcid_get_vf3, funcid_get_int_array_length, funcid_set_int_array, funcid_get_int_array,
    funcid_get_float_array_length, funcid_set_float_array, funcid_get_float_array,
    funcid_get_vf3_array_length, funcid_set_vf3_array, funcid_get_vf3_array
};

#endif
