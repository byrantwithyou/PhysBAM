//#####################################################################
// Copyright 2004-2008, Geoffrey Irving, Frank Losasso, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG_PRINTF.h>
#include <Core/Log/SCOPE.h>
#include <cassert>
#include <ostream>
#include <sstream>
namespace PhysBAM{
namespace LOG{
//#####################################################################
// Constructor
//#####################################################################
SCOPE::
SCOPE()
    :active(false)
{
}
//#####################################################################
// Constructor
//#####################################################################
SCOPE::
SCOPE(const std::string& scope_identifier)
    :active(true)
{
    LOG_CLASS::Push_Scope(scope_identifier,scope_identifier);
}
//#####################################################################
// Constructor
//#####################################################################
SCOPE::
SCOPE(const std::string& scope_identifier,const std::string& scope_name)
    :active(true)
{
    LOG_CLASS::Push_Scope(scope_identifier,scope_name);
}
//#####################################################################
// Destructor
//#####################################################################
SCOPE::
~SCOPE()
{
    if(active) LOG_CLASS::Pop_Scope();
}
//#####################################################################
// Function Push
//#####################################################################
void SCOPE::
Push(const std::string& scope_identifier,const std::string& scope_name)
{
    assert(!active);
    active=true;
    LOG_CLASS::Push_Scope(scope_identifier,scope_name);
}
//#####################################################################
// Function Pop
//#####################################################################
void SCOPE::
Pop()
{
    assert(active);
    active=false;
    LOG_CLASS::Pop_Scope();
}
//##################################################################### 
}
}

