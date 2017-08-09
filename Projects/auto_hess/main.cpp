//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include "helper.h"

using namespace PhysBAM;

void Test_S();
void Test_V();
void Test_M();
void Test_MV();
void Test_MM();
void Test_MM_mult();

int main(int argc, char* argv[])
{
    Test_S();
    Test_V();
    Test_M();
    Test_MV();
    Test_MM();
    Test_MM_mult();
    
    return 0;
}

