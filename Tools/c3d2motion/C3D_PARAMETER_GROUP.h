//##########################################################################################
// Class C3D_PARAMETER_GROUP
// A group in the param table, holds a list of values
//##########################################################################################
#ifndef _C3D_PARAMETER_GROUP_h
#define C3D_PARAMETER_GROUP_h
#include "C3D_PARAMETER_VALUE.h"

class C3D_PARAMETER_GROUP:public C3D_PARAMETER_ENTRY{
    std::map<std::string,C3D_PARAMETER_VALUE*> values;
public:
    C3D_PARAMETER_GROUP()
        :C3D_PARAMETER_ENTRY() 
    {}
    
    virtual ~C3D_PARAMETER_GROUP()
    {std::map<std::string,C3D_PARAMETER_VALUE*>::iterator i;for(i=values.begin();i!=values.end();++i){delete i->second;}}
    
    void Display_Data()
    {std::map<std::string,C3D_PARAMETER_VALUE*>::iterator i;for (i=values.begin();i!=values.end();++i) i->second->Display_Data();}
    
    void Add_Member(const char* name,C3D_PARAMETER_VALUE* entry)
    {values[std::string(name)]=entry;}
    
    C3D_PARAMETER_VALUE* Get(const char *name)
    {return values[name];}
    
//##########################################################################################
};

#endif
