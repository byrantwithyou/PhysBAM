//##########################################################################################
// Class  C3D_PARAMETER_ENTRY
// An entry in the parameter table (a value or group,see below)
//##########################################################################################
#ifndef _C3D_PARAMETER_ENTRY_h
#define C3D_PARAMETER_ENTRY_h
#include <string>
#include <stdlib.h>
#include <string.h>

class C3D_PARAMETER_ENTRY{
    std::string name,description;
public:
    C3D_PARAMETER_ENTRY(){}
    virtual ~C3D_PARAMETER_ENTRY(){}
    
    virtual void Set_Info(const char* name_input,const char* description_input)
    {name=name_input;description=description_input;}
    
    const char* Get_Name(){return name.c_str();}
    const char* Get_Description(){return description.c_str();}
};

#endif
