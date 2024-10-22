//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// File INTERRUPTS
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Utilities/INTERRUPTS.h>
namespace PhysBAM{
//#####################################################################
namespace {
ARRAY<void(*)()> interrupt_checkers;
}
//#####################################################################
// Function Check_For_Interrupts
//#####################################################################
void Check_For_Interrupts()
{
    for(int i=0;i<interrupt_checkers.Size();i++)
        interrupt_checkers(i)();
}
//#####################################################################
// Function Check_For_Interrupts
//#####################################################################
void Add_Interrupt_Checker(void (*checker)())
{
    interrupt_checkers.Append(checker);
}
//#####################################################################
}
