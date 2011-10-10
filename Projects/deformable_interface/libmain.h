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
int get_array_length(const physbam_base * obj, int id);
void set_int_array(physbam_base * obj, int id, const int * x, int length, int start = 0);
void get_int_array(const physbam_base * obj, int id, int * x, int length, int start = 0);
void set_float_array(physbam_base * obj, int id, const float * x, int length, int start = 0);
void get_float_array(const physbam_base * obj, int id, float * x, int length, int start = 0);
}

#endif
