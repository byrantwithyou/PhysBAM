//#####################################################################
// Copyright 2006, Geoffrey Irving, Eftychios Sifakis, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Core/Read_Write/TYPED_STREAM.h>
#include <iostream>
#include <memory>
using namespace PhysBAM;
//#####################################################################
// Function Set
//#####################################################################
void FILE_ISTREAM::
Set(std::istream* stream_input)
{
    stream=stream_input;
    Read_Primitive(*stream,type.use_doubles);
}
//#####################################################################
// Destructor
//#####################################################################
FILE_ISTREAM::
~FILE_ISTREAM()
{
    delete stream;
}
//#####################################################################
// Function Set
//#####################################################################
void FILE_OSTREAM::
Set(std::ostream* stream_input,STREAM_TYPE type_input)
{
    stream=stream_input;
    type=type_input;
    Write_Primitive(*stream,type.use_doubles);
}
//#####################################################################
// Destructor
//#####################################################################
FILE_OSTREAM::
~FILE_OSTREAM()
{
    delete stream;
}
