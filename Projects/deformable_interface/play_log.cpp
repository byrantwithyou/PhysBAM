#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include "DATA_EXCHANGE_CONVERSION.h"
#include "libmain.h"
#include "request_logging.h"
#include <dlfcn.h>

struct physbam_simulation;
struct physbam_object;
struct physbam_force;

struct simulation_info
{
    physbam_simulation* sim;
    std::map<int,physbam_base*> objects;
};

int main(int argc,char* argv[])
{
    using namespace data_exchange;
    register_ids();
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
    int (*get_id)(physbam_base * obj, const char * attribute) = (int (*)(physbam_base * obj, const char * attribute)) dlsym(handle, "get_id");
    void (*set_int)(physbam_base * obj, int id, int x) = (void (*)(physbam_base * obj, int id, int x)) dlsym(handle, "set_int");
    int (*get_int)(const physbam_base * obj, int id) = (int (*)(const physbam_base * obj, int id)) dlsym(handle, "get_int");
    void (*set_float)(physbam_base * obj, int id, float x) = (void (*)(physbam_base * obj, int id, float x)) dlsym(handle, "set_float");
    float (*get_float)(const physbam_base * obj, int id) = (float (*)(const physbam_base * obj, int id)) dlsym(handle, "get_float");

    const char* log = argv[1];
    if(!log){
        printf("Usage: %s <log>\n", argv[0]);
        exit(0);
    }

    std::ifstream in(log);
    std::map<int,simulation_info> simulations;

    while(in)
    {
        int func_id = -1;
        int sim_id = -1;
        PhysBAM::Read_Binary<float>(in,func_id,sim_id);
        simulation_info& si = simulations[sim_id];
        switch(func_id)
        {
            case funcid_create_simulation:
                si.sim = create_simulation();
                break;

            case funcid_destroy_simulation:
                destroy_simulation(si.sim);
                simulations.erase(sim_id);
                break;

            case funcid_add_object:{
                data_exchange::simulation_object* so = 0;
                Read(in,so);
                int so_id = -1;
                PhysBAM::Read_Binary<float>(in,so_id);
                si.objects[so_id] = add_object(si.sim,so);
                delete so;
                break;}

            case funcid_add_force:{
                data_exchange::force* f = 0;
                Read(in,f);
                int f_id = -1;
                PhysBAM::Read_Binary<float>(in,f_id);
                si.objects[f_id] = add_force(si.sim,f);
                delete f;
                break;}

            case funcid_apply_force_to_object:{
                int oi = -1, fi = -1;
                PhysBAM::Read_Binary<float>(in,oi,fi);
                apply_force_to_object((physbam_object*)si.objects[oi],(physbam_force*)si.objects[fi]);
                break;}

            case funcid_simulate_frame:
                simulate_frame(si.sim);
                break;

            case funcid_get_id:{
                int oi = -1;
                std::string attribute;
                PhysBAM::Read_Binary<float>(in,oi,attribute);
                get_id(si.objects[oi], attribute.c_str());
                break;}
 
            case funcid_set_int:{
                int oi = -1, id = -1, x;
                PhysBAM::Read_Binary<float>(in,oi,id,x);
                set_int(si.objects[oi], id, x);
                break;}

            case funcid_get_int:{
                int oi = -1, id = -1;
                PhysBAM::Read_Binary<float>(in,oi,id);
                get_int(si.objects[oi], id);
                break;}

            case funcid_set_float:{
                int oi = -1, id = -1;
                float x;
                PhysBAM::Read_Binary<float>(in,oi,id,x);
                set_float(si.objects[oi], id, x);
                break;}

            case funcid_get_float:{
                int oi = -1, id = -1;
                PhysBAM::Read_Binary<float>(in,oi,id);
                get_float(si.objects[oi], id);
                break;}

            case funcid_set_vf3:{
                int oi = -1, id = -1;
                PhysBAM::VECTOR<float,3> x;
                PhysBAM::Read_Binary<float>(in,oi,id,x);
                data_exchange::vf3 v;
                From_Pb(v, x);
                set_vf3(si.objects[oi], id, v);
                break;}

            case funcid_get_vf3:{
                int oi = -1, id = -1;
                PhysBAM::Read_Binary<float>(in,oi,id);
                get_vf3(si.objects[oi], id);
                break;}

            case funcid_get_int_array_length:{
                int oi = -1, id = -1;
                PhysBAM::Read_Binary<float>(in,oi,id);
                get_int_array_length(si.objects[oi], id);
                break;}

            case funcid_set_int_array:{
                int oi = -1, id = -1, start = -1;
                ARRAY<int> array;
                PhysBAM::Read_Binary<float>(in,oi,id,array,start);
                set_int_array(si.objects[oi], id, array.Get_Array_Pointer(), array.m, start);
                break;}

            case funcid_get_int_array:{
                int oi = -1, id = -1, length = -1, start = -1;
                PhysBAM::Read_Binary<float>(in,oi,id,length,start);
                ARRAY<int> array(length);
                get_int_array(si.objects[oi], id, array.Get_Array_Pointer(), length, start);
                break;}

            case funcid_get_float_array_length:{
                int oi = -1, id = -1;
                PhysBAM::Read_Binary<float>(in,oi,id);
                get_float_array_length(si.objects[oi], id);
                break;}

            case funcid_set_float_array:{
                int oi = -1, id = -1, start = -1;
                ARRAY<float> array;
                PhysBAM::Read_Binary<float>(in,oi,id,array,start);
                set_float_array(si.objects[oi], id, array.Get_Array_Pointer(), array.m, start);
                break;}

            case funcid_get_float_array:{
                int oi = -1, id = -1, length = -1, start = -1;
                PhysBAM::Read_Binary<float>(in,oi,id,length,start);
                ARRAY<float> array(length);
                get_float_array(si.objects[oi], id, array.Get_Array_Pointer(), length, start);
                break;}

            case funcid_get_vf3_array_length:{
                int oi = -1, id = -1;
                PhysBAM::Read_Binary<float>(in,oi,id);
                get_vf3_array_length(si.objects[oi], id);
                break;}

            case funcid_set_vf3_array:{
                int oi = -1, id = -1, start = -1;
                ARRAY<VECTOR<float,3> > array;
                PhysBAM::Read_Binary<float>(in,oi,id,array,start);
                set_vf3_array(si.objects[oi], id, (data_exchange::vf3*) array.Get_Array_Pointer(), array.m, start);
                break;}

            case funcid_get_vf3_array:{
                int oi = -1, id = -1, length = -1, start = -1;
                PhysBAM::Read_Binary<float>(in,oi,id,length,start);
                ARRAY<VECTOR<float,3> > array(length);
                get_vf3_array(si.objects[oi], id, (data_exchange::vf3*) array.Get_Array_Pointer(), length, start);
                break;}

            default: printf("Unrecognized function id: %i\n", func_id); exit(1);
        }
    }

    

    return 0;
}

