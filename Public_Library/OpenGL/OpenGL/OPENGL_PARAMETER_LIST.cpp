//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// OpenGL-related functions for PARAMETER_LIST
//##################################################################### 
#include <Core/Read_Write/STRING_UTILITIES.h>
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <sstream>
namespace PhysBAM{
//#####################################################################
// Parsing for OPENGL_COLOR
//#####################################################################
template<> void PARAMETER_LIST::
Value_To_String(const OPENGL_COLOR &value,std::string &value_string)
{
    std::ostringstream value_string_stream;
    value_string_stream << value;
    value_string = value_string_stream.str();
}
template<> bool PARAMETER_LIST::
String_To_Value(const std::string &value_string,OPENGL_COLOR &value,bool commandline_style)
{
    std::istringstream value_string_stream(value_string);
    std::string str;
    value_string_stream >> str;
    if(Is_Number(str)){
        String_To_Value(str,value.rgba[0]);
        value_string_stream >> value.rgba[1] >> value.rgba[2];}
    else{ // case insensitive
        if(!Compare_Strings(str,"Red",false)) value=OPENGL_COLOR::Red();
        else if(!Compare_Strings(str,"Green",false)) value=OPENGL_COLOR::Green();
        else if(!Compare_Strings(str,"Blue",false)) value=OPENGL_COLOR::Blue();
        else if(!Compare_Strings(str,"Yellow",false)) value=OPENGL_COLOR::Yellow();
        else if(!Compare_Strings(str,"Cyan",false)) value=OPENGL_COLOR::Cyan();
        else if(!Compare_Strings(str,"Magenta",false)) value=OPENGL_COLOR::Magenta();
        else{LOG::cerr<<"Unrecognized color '"<<str<<"'"<<std::endl;return false;}}
    bool success = !value_string_stream.fail();
    if(!(value_string_stream >> value.rgba[3])) {
        value.rgba[3] = 1;
    }
    return success;
}
//#####################################################################
};
