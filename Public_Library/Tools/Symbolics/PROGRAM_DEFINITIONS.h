//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
// Class PROGRAM_DEFINITIONS
//#####################################################################
#ifndef __PROGRAM_DEFINITIONS__
#define __PROGRAM_DEFINITIONS__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Symbolics/CODE_BLOCK.h>
#include <Tools/Symbolics/CODE_BLOCK_NODE.h>
#include <Tools/Symbolics/INSTRUCTION.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
#include <map>
#include <string>

namespace PhysBAM{

struct PROGRAM_PARSE_NODE
{
    int type;
    int val;
    PROGRAM_PARSE_NODE *a,*b;

    PROGRAM_PARSE_NODE(int type,int val,PROGRAM_PARSE_NODE *a,PROGRAM_PARSE_NODE *b)
        :type(type),val(val),a(a),b(b)
    {}

    ~PROGRAM_PARSE_NODE(){delete a;delete b;}
};

extern ARRAY<std::string> parse_identifiers;
extern ARRAY<double> parse_constants;
extern PROGRAM_PARSE_NODE* parse_root;

}
#endif
