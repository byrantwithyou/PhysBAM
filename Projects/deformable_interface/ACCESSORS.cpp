#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <string>
#include "DEFORMABLE_EXAMPLE.h"
using namespace PhysBAM;
ARRAY<HASHTABLE<std::string,int> > attribute_name_lookup;

#define CREATE_LOOKUP_DATABASE(attribute_type) \
    struct attribute_accessor_##attribute_type \
    { \
        attribute_type (*get)(const BASE_WRAPPER& w); \
        void (*set)(BASE_WRAPPER& w, attribute_type x); \
        attribute_accessor_##attribute_type() {} \
        attribute_accessor_##attribute_type(attribute_type (*g)(const BASE_WRAPPER& w),void (*s)(BASE_WRAPPER& w, attribute_type x)): get(g), set(s) {} \
    }; \
    ARRAY<ARRAY<attribute_accessor_##attribute_type> > attribute_accessor_lookup_##attribute_type

CREATE_LOOKUP_DATABASE(int);
CREATE_LOOKUP_DATABASE(float);

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


void PhysBAM::register_accessors(int last_id)
{
    attribute_name_lookup.Resize(last_id);
    attribute_accessor_lookup_int.Resize(last_id);
    attribute_accessor_lookup_float.Resize(last_id);

    REGISTER_ACCESSOR_TYPE(GRAVITY_WRAPPER,magnitude,float,{return w.magnitude;},{
            w.magnitude=x;
            for(int i=1;i<=w.force_instances.m;i++)
                w.force_instances(i)->gravity=x;
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

