#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <string>
#include "DATA_EXCHANGE_CONVERSION.h"
#include "DEFORMABLE_EXAMPLE.h"
using namespace PhysBAM;
typedef VECTOR<float,3> TV;
ARRAY<HASHTABLE<std::string,int> > attribute_name_lookup;

#define CREATE_LOOKUP_DATABASE(attribute_type) \
    struct attribute_accessor_##attribute_type \
    { \
        attribute_type (*get)(const BASE_WRAPPER& w); \
        void (*set)(BASE_WRAPPER& w, attribute_type x); \
        attribute_accessor_##attribute_type(): get(0), set(0) {}        \
        attribute_accessor_##attribute_type(attribute_type (*g)(const BASE_WRAPPER&),void (*s)(BASE_WRAPPER&, attribute_type)): get(g), set(s) {} \
    }; \
    ARRAY<ARRAY<attribute_accessor_##attribute_type> > attribute_accessor_lookup_##attribute_type

#define CREATE_LOOKUP_ARRAY_DATABASE(attribute_type) \
    struct attribute_accessor_array_##attribute_type \
    { \
        typedef int (*Z)(const BASE_WRAPPER& w); \
        typedef void (*G)(const BASE_WRAPPER& w,ARRAY_VIEW<attribute_type>&,int offset); \
        typedef void (*S)(BASE_WRAPPER& w,const ARRAY_VIEW<const attribute_type>&,int offset); \
        Z size; \
        G get; \
        S set; \
        attribute_accessor_array_##attribute_type(): size(0), get(0), set(0) {} \
        attribute_accessor_array_##attribute_type(Z z,G g,S s): size(z), get(g), set(s) {} \
    }; \
    ARRAY<ARRAY<attribute_accessor_array_##attribute_type> > attribute_accessor_lookup_array_##attribute_type

CREATE_LOOKUP_DATABASE(int);
CREATE_LOOKUP_DATABASE(float);
CREATE_LOOKUP_DATABASE(TV);
CREATE_LOOKUP_ARRAY_DATABASE(int);
CREATE_LOOKUP_ARRAY_DATABASE(float);
CREATE_LOOKUP_ARRAY_DATABASE(TV);

#define REGISTER_ACCESSOR_TYPE(class_name,attribute_name,attribute_type,get_body,set_body) \
    struct class_name##_##attribute_name \
    { \
        static attribute_type get(const class_name& w) get_body \
        static void set(class_name& w, attribute_type x) set_body \
    }; \
    int class_name##_##attribute_name##_id = attribute_accessor_lookup_##attribute_type(class_name::fixed_id()).Append(attribute_accessor_##attribute_type( \
        (attribute_type (*)(const BASE_WRAPPER&))class_name##_##attribute_name::get, \
        (void (*)(BASE_WRAPPER&,attribute_type))class_name##_##attribute_name::set)); \
    attribute_name_lookup(class_name::fixed_id()).Set(#attribute_name,class_name##_##attribute_name##_id);

#define REGISTER_ACCESSOR_ARRAY_TYPE(class_name,attribute_name,attribute_type,size_body,get_body,set_body) \
    struct class_name##_##attribute_name \
    { \
        static int size(const class_name& w) size_body \
        static void get(const class_name& w,ARRAY_VIEW<attribute_type>& array,int offset) get_body \
        static void set(class_name& w,const ARRAY_VIEW<const attribute_type>& array,int offset) set_body \
    }; \
    int class_name##_##attribute_name##_id = attribute_accessor_lookup_array_##attribute_type(class_name::fixed_id()).Append(attribute_accessor_array_##attribute_type( \
        (int (*)(const BASE_WRAPPER&))class_name##_##attribute_name::size, \
        (void (*)(const BASE_WRAPPER&,ARRAY_VIEW<attribute_type>&,int))class_name##_##attribute_name::get, \
        (void (*)(BASE_WRAPPER&,const ARRAY_VIEW<const attribute_type>&,int))class_name##_##attribute_name::set)); \
    attribute_name_lookup(class_name::fixed_id()).Set(#attribute_name,class_name##_##attribute_name##_id);

//int get_array_length(const physbam_base * obj, int id);
//void set_float_array(physbam_base * obj, int id, const float * x, int length, int start = 0);
//void get_float_array(const physbam_base * obj, int id, float * x, int length, int start = 0);

void PhysBAM::register_accessors(int last_id)
{
    attribute_name_lookup.Resize(last_id);
    attribute_accessor_lookup_int.Resize(last_id);
    attribute_accessor_lookup_float.Resize(last_id);
    attribute_accessor_lookup_TV.Resize(last_id);
    attribute_accessor_lookup_array_int.Resize(last_id);
    attribute_accessor_lookup_array_float.Resize(last_id);
    attribute_accessor_lookup_array_TV.Resize(last_id);

    // GRAVITY_WRAPPER

    REGISTER_ACCESSOR_TYPE(GRAVITY_WRAPPER,magnitude,float,{return w.magnitude;},{
            w.magnitude=x;
            for(int i=1;i<=w.force_instances.m;i++)
                w.force_instances(i)->gravity=x;
        });

    REGISTER_ACCESSOR_TYPE(GRAVITY_WRAPPER,direction,TV,{return w.direction;},{
            w.direction=x;
            for(int i=1;i<=w.force_instances.m;i++)
                w.force_instances(i)->downward_direction=x;
        });

    // VOLUMETRIC_FORCE_WRAPPER

    REGISTER_ACCESSOR_TYPE(VOLUMETRIC_FORCE_WRAPPER,stiffness,float,{return w.stiffness;},{w.stiffness=x;w.Propagate_Parameters();});
    REGISTER_ACCESSOR_TYPE(VOLUMETRIC_FORCE_WRAPPER,poissons_ratio,float,{return w.poissons_ratio;},{w.poissons_ratio=x;w.Propagate_Parameters();});
    REGISTER_ACCESSOR_TYPE(VOLUMETRIC_FORCE_WRAPPER,damping,float,{return w.damping;},{w.damping=x;w.Propagate_Parameters();});

    // DEFORMABLE_BODY_WRAPPER

    REGISTER_ACCESSOR_ARRAY_TYPE(DEFORMABLE_BODY_WRAPPER,position,TV,{return w.particle_map.m;},
        {
            for(int i=1;i<=array.m;i++) array(i)=w.de.solid_body_collection.deformable_body_collection.particles.X(w.particle_map(i+offset));
        },
        {
            for(int i=1;i<=array.m;i++) w.de.solid_body_collection.deformable_body_collection.particles.X(w.particle_map(i+offset))=array(i);
        });

    REGISTER_ACCESSOR_ARRAY_TYPE(DEFORMABLE_BODY_WRAPPER,velocity,TV,{return w.particle_map.m;},
        {
            for(int i=1;i<=array.m;i++) array(i)=w.de.solid_body_collection.deformable_body_collection.particles.V(w.particle_map(i+offset));
        },
        {
            for(int i=1;i<=array.m;i++) w.de.solid_body_collection.deformable_body_collection.particles.V(w.particle_map(i+offset))=array(i);
        });


}

int get_id(physbam_base * obj, const char * attribute)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    int * id = attribute_name_lookup(w->id).Get_Pointer(attribute);
    if(id) return *id;
    return -1;
}

void set_int(physbam_base * obj, int id, int x)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    attribute_accessor_lookup_int(w->id)(id).set(*w,x);
}

int get_int(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    return attribute_accessor_lookup_int(w->id)(id).get(*w);
}

void set_float(physbam_base * obj, int id, float x)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    attribute_accessor_lookup_float(w->id)(id).set(*w,x);
}

float get_float(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    return attribute_accessor_lookup_float(w->id)(id).get(*w);
}

void set_vf3(physbam_base * obj, int id, data_exchange::vf3 x)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    attribute_accessor_lookup_TV(w->id)(id).set(*w,To_Pb(x));
}

data_exchange::vf3 get_vf3(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    data_exchange::vf3 ret;
    From_Pb(ret,attribute_accessor_lookup_TV(w->id)(id).get(*w));
    return ret;
}

int get_int_array_length(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    return attribute_accessor_lookup_array_int(w->id)(id).size(*w);
}

void set_int_array(physbam_base * obj, int id, const int * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<const int> array(length,x);
    attribute_accessor_lookup_array_int(w->id)(id).set(*w,array,start);
}

void get_int_array(const physbam_base * obj, int id, int * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<int> array(length,x);
    return attribute_accessor_lookup_array_int(w->id)(id).get(*w,array,start);
}

int get_float_array_length(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    return attribute_accessor_lookup_array_float(w->id)(id).size(*w);
}

void set_float_array(physbam_base * obj, int id, const float * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<const float> array(length,x);
    attribute_accessor_lookup_array_float(w->id)(id).set(*w,array,start);
}

void get_float_array(const physbam_base * obj, int id, float * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<float> array(length,x);
    return attribute_accessor_lookup_array_float(w->id)(id).get(*w,array,start);
}

int get_vf3_array_length(const physbam_base * obj, int id)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    return attribute_accessor_lookup_array_TV(w->id)(id).size(*w);
}

void set_vf3_array(physbam_base * obj, int id, const data_exchange::vf3 * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<const TV> array(length,(TV*) x);
    attribute_accessor_lookup_array_TV(w->id)(id).set(*w,array,start);
}

void get_vf3_array(const physbam_base * obj, int id, data_exchange::vf3 * x, int length, int start)
{
    BASE_WRAPPER * w = (BASE_WRAPPER *) obj;
    ARRAY_VIEW<TV> array(length,(TV*) x);
    return attribute_accessor_lookup_array_TV(w->id)(id).get(*w,array,start);
}

