//##########################################################################################
// Class  C3D_PARSE_ERROR
// what happens when bad things go down
//##########################################################################################
#ifndef _C3D_PARSE_ERROR_h
#define _C3D_PARSE_ERROR_h
#include <string>

class C3D_PARSE_ERROR{
public:
    std::string error_string;
    C3D_PARSE_ERROR(std::string error_string_input)
        :error_string(error_string_input)
    {}
};

#endif
