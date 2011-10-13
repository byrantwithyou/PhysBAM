#ifndef __ACCESSORS__
#define __ACCESSORS__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <string>
#include "DEFORMABLE_EXAMPLE.h"
using namespace PhysBAM;
extern ARRAY<HASHTABLE<std::string,int> > attribute_name_lookup;
typedef VECTOR<float,3> TV;

#define CREATE_LOOKUP_DATABASE(attribute_type) \
    struct attribute_accessor_##attribute_type \
    { \
        attribute_type (*get)(const BASE_WRAPPER& w); \
        void (*set)(BASE_WRAPPER& w, attribute_type x); \
        attribute_accessor_##attribute_type(): get(0), set(0) {}        \
        attribute_accessor_##attribute_type(attribute_type (*g)(const BASE_WRAPPER&),void (*s)(BASE_WRAPPER&, attribute_type)): get(g), set(s) {} \
    }; \
    extern ARRAY<ARRAY<attribute_accessor_##attribute_type> > attribute_accessor_lookup_##attribute_type

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
    extern ARRAY<ARRAY<attribute_accessor_array_##attribute_type> > attribute_accessor_lookup_array_##attribute_type


#define CREATE_LOOKUP_DATABASE(attribute_type) \
    struct attribute_accessor_##attribute_type \
    { \
        attribute_type (*get)(const BASE_WRAPPER& w); \
        void (*set)(BASE_WRAPPER& w, attribute_type x); \
        attribute_accessor_##attribute_type(): get(0), set(0) {}        \
        attribute_accessor_##attribute_type(attribute_type (*g)(const BASE_WRAPPER&),void (*s)(BASE_WRAPPER&, attribute_type)): get(g), set(s) {} \
    }; \
    extern ARRAY<ARRAY<attribute_accessor_##attribute_type> > attribute_accessor_lookup_##attribute_type

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
    extern ARRAY<ARRAY<attribute_accessor_array_##attribute_type> > attribute_accessor_lookup_array_##attribute_type

CREATE_LOOKUP_DATABASE(int);
CREATE_LOOKUP_DATABASE(float);
CREATE_LOOKUP_DATABASE(TV);
CREATE_LOOKUP_ARRAY_DATABASE(int);
CREATE_LOOKUP_ARRAY_DATABASE(float);
CREATE_LOOKUP_ARRAY_DATABASE(TV);

#endif
