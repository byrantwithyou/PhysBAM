//#####################################################################
// Copyright 2004-2008, Geoffrey Irving, Frank Losasso, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCOPE
//##################################################################### 
#ifndef __SCOPE__
#define __SCOPE__

#include <Tools/Log/LOG_PRINTF.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <cassert>
#include <ostream>
#include <sstream>
namespace PhysBAM{
namespace LOG{
class SCOPE:private NONCOPYABLE
{
    bool active;
public:
    SCOPE();
    SCOPE(const std::string& scope_identifier);
    SCOPE(const std::string& scope_identifier,const std::string& scope_name);

    template<class T1,class ...Args>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1,Args&& ...args)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,LOG::sprintf(format.c_str(),d1,args...));
    }

    ~SCOPE();
    
    void Push(const std::string& scope_identifier,const std::string& scope_name);

    template<class T1>
    void Push(const std::string& scope_identifier,const std::string& format,const T1& d1)
    {Push(scope_identifier,LOG::sprintf(format,d1));}

    void Pop();
//##################################################################### 
};
}
}
#endif
