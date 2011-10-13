#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <string>
#include "ACCESSORS.h"
using namespace PhysBAM;

ARRAY<HASHTABLE<std::string,int> > attribute_name_lookup;
ARRAY<ARRAY<attribute_accessor_int> > attribute_accessor_lookup_int;
ARRAY<ARRAY<attribute_accessor_float> > attribute_accessor_lookup_float;
ARRAY<ARRAY<attribute_accessor_TV> > attribute_accessor_lookup_TV;
ARRAY<ARRAY<attribute_accessor_array_int> > attribute_accessor_lookup_array_int;
ARRAY<ARRAY<attribute_accessor_array_float> > attribute_accessor_lookup_array_float;
ARRAY<ARRAY<attribute_accessor_array_TV> > attribute_accessor_lookup_array_TV;

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
