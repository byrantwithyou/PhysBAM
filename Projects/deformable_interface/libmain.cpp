//#####################################################################
// Copyright 2011, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include "ACCESSORS.h"
#include "DATA_EXCHANGE_CONVERSION.h"
#include "DEFORMABLE_EXAMPLE.h"
#include "libmain.h"
#include "request_logging.h"
#include <stdio.h>
using namespace PhysBAM;
using namespace data_exchange;

void PhysBAM::Register_Wrapper_Ids();

physbam_simulation * create_simulation()
{
    data_exchange::register_ids();
    Register_Wrapper_Ids();
    DEFORMABLE_EXAMPLE<float> * de = new DEFORMABLE_EXAMPLE<float>(STREAM_TYPE(0.f));
    de->Initialize_Simulation();
    if(log_info)
    {
        Write_Binary<float>(log_info->out,funcid_create_simulation,log_info->next_simulation_id);
        log_info->simulation_map[(physbam_simulation*) de].id = log_info->next_simulation_id++;
    }
    return (physbam_simulation*) de;
}

bool destroy_simulation(physbam_simulation * sim)
{
    DEFORMABLE_EXAMPLE<float> * de = (DEFORMABLE_EXAMPLE<float>*) sim;
    delete de;
    if(log_info)
    {
        Write_Binary<float>(log_info->out,funcid_destroy_simulation,log_info->simulation_map[sim].id);
        log_info->simulation_map.erase(sim);
    }
    return true;
}

physbam_object * add_object(physbam_simulation * sim, const simulation_object* so)
{
    DEFORMABLE_EXAMPLE<float> * de = (DEFORMABLE_EXAMPLE<float>*) sim;
    OBJECT_WRAPPER* wrap=de->Add_Simulation_Object(*so);
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[sim];
        Write_Binary<float>(log_info->out,funcid_add_object,si.id);
        Write(log_info->out,so);
        Write_Binary<float>(log_info->out,si.next_object);
        si.object_map[(physbam_object*) wrap] = si.next_object++;
    }
    return (physbam_object*) wrap;
}

physbam_force * add_force(physbam_simulation * sim, const force* f)
{
    DEFORMABLE_EXAMPLE<float> * de = (DEFORMABLE_EXAMPLE<float>*) sim;
    FORCE_WRAPPER* wrap=de->Add_Force(*f);
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[sim];
        Write_Binary<float>(log_info->out,funcid_add_force,si.id);
        Write(log_info->out,f);
        Write_Binary<float>(log_info->out,si.next_object);
        si.object_map[(physbam_force*) wrap] = si.next_object++;
    }
    return (physbam_force*) wrap;
}

bool apply_force_to_object(physbam_object * obj, physbam_force* f)
{
    OBJECT_WRAPPER * ow = (OBJECT_WRAPPER *) obj;
    FORCE_WRAPPER * fw = (FORCE_WRAPPER *) f;
    ow->de.new_forces_relations.Append(PAIR<OBJECT_WRAPPER*,FORCE_WRAPPER*>(ow,fw));
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &ow->de];
        Write_Binary<float>(log_info->out,funcid_apply_force_to_object,si.id,si.object_map[obj],si.object_map[f]);
    }
    return true;
}

bool simulate_frame(physbam_simulation * sim)
{
    DEFORMABLE_EXAMPLE<float> * de = (DEFORMABLE_EXAMPLE<float>*) sim;
    if(log_info)
    {
        Write_Binary<float>(log_info->out,funcid_simulate_frame,log_info->simulation_map[sim].id);
    }
    de->Simulate_Frame();
    return true;
}

int get_id(physbam_base * obj, const char * attribute)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_id,si.id,si.object_map[obj],attribute);
    }
    int * id = attribute_name_lookup(w->id).Get_Pointer(attribute);
    if(id) return *id;
    return -1;
}

void set_int(physbam_base * obj, int id, int x)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_set_int,si.id,si.object_map[obj],id,x);
    }
    attribute_accessor_lookup_int(w->id)(id).set(*w,x);
}

int get_int(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_int,si.id,si.object_map[obj],id);
    }
    return attribute_accessor_lookup_int(w->id)(id).get(*w);
}

void set_float(physbam_base * obj, int id, float x)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_set_float,si.id,si.object_map[obj],id,x);
    }
    attribute_accessor_lookup_float(w->id)(id).set(*w,x);
}

float get_float(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_float,si.id,si.object_map[obj],id);
    }
    return attribute_accessor_lookup_float(w->id)(id).get(*w);
}

void set_vf3(physbam_base * obj, int id, data_exchange::vf3 x)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_set_vf3,si.id,si.object_map[obj],id);
        Write(log_info->out,x);
    }
    attribute_accessor_lookup_TV(w->id)(id).set(*w,To_Pb(x));
}

data_exchange::vf3 get_vf3(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    data_exchange::vf3 ret;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_vf3,si.id,si.object_map[obj],id);
    }
    From_Pb(ret,attribute_accessor_lookup_TV(w->id)(id).get(*w));
    return ret;
}

int get_int_array_length(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_int_array_length,si.id,si.object_map[obj],id);
    }
    return attribute_accessor_lookup_array_int(w->id)(id).size(*w);
}

void set_int_array(physbam_base * obj, int id, const int * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<const int> array(length,x);
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_set_int_array,si.id,si.object_map[obj],id,array,start);
    }
    attribute_accessor_lookup_array_int(w->id)(id).set(*w,array,start);
}

void get_int_array(const physbam_base * obj, int id, int * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<int> array(length,x);
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_int_array,si.id,si.object_map[obj],id,length,start);
    }
    return attribute_accessor_lookup_array_int(w->id)(id).get(*w,array,start);
}

int get_float_array_length(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_float_array_length,si.id,si.object_map[obj],id);
    }
    return attribute_accessor_lookup_array_float(w->id)(id).size(*w);
}

void set_float_array(physbam_base * obj, int id, const float * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<const float> array(length,x);
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_set_float_array,si.id,si.object_map[obj],id,array,start);
    }
    attribute_accessor_lookup_array_float(w->id)(id).set(*w,array,start);
}

void get_float_array(const physbam_base * obj, int id, float * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<float> array(length,x);
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_float_array,si.id,si.object_map[obj],id,length,start);
    }
    return attribute_accessor_lookup_array_float(w->id)(id).get(*w,array,start);
}

int get_vf3_array_length(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_vf3_array_length,si.id,si.object_map[obj],id);
    }
    return attribute_accessor_lookup_array_TV(w->id)(id).size(*w);
}

void set_vf3_array(physbam_base * obj, int id, const data_exchange::vf3 * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<const TV> array(length,(TV*) x);
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_set_vf3_array,si.id,si.object_map[obj],id,array,start);
    }
    attribute_accessor_lookup_array_TV(w->id)(id).set(*w,array,start);
}

void get_vf3_array(const physbam_base * obj, int id, data_exchange::vf3 * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<TV> array(length,(TV*) x);
    if(log_info)
    {
        request_log_simulation_info& si = log_info->simulation_map[(physbam_simulation*) &w->de];
        Write_Binary<float>(log_info->out,funcid_get_vf3_array,si.id,si.object_map[obj],id,length,start);
    }
    return attribute_accessor_lookup_array_TV(w->id)(id).get(*w,array,start);
}


