//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Solids/BW_DRIVER.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "2d_Tests/STANDARD_TESTS_2D.h"
#include "Baraff_Witkin_Tests/BARAFF_WITKIN_TESTS.h"
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Toy_Problems/TOY_PROBLEMS.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    if(PARSE_ARGS::Find_And_Remove("-2d",argc,argv)){
        typedef VECTOR<T,2> TV;
        EXAMPLE<TV>* example;
        DRIVER<TV>* driver;

        example=new STANDARD_TESTS_2D<T>(stream_type);
        example->Parse(argc,argv);
        driver=new SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<VECTOR<T,2> > >(*dynamic_cast<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,2> > >*>(example));
        driver->Execute_Main_Program();
        
        delete driver;delete example;}
    else{
        typedef VECTOR<T,3> TV;
        EXAMPLE<TV>* example;
        DRIVER<TV>* driver;

        if(PARSE_ARGS::Find_And_Remove("-bw",argc,argv)){
            example=new BARAFF_WITKIN_TESTS<T>(stream_type);
            example->Parse(argc,argv);
            driver=new BW_DRIVER<TV>(*dynamic_cast<SOLIDS_EXAMPLE<TV>*>(example));}
        else if(PARSE_ARGS::Find_And_Remove("-toy",argc,argv)){
            example=new TOY_PROBLEMS<T>(stream_type);
            example->Parse(argc,argv);
            driver=new SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> >(*dynamic_cast<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >*>(example));}
        else{
            example=new STANDARD_TESTS<T>(stream_type);
            example->Parse(argc,argv);
            driver=new SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> >(*dynamic_cast<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >*>(example));}
        driver->Execute_Main_Program();
        
        delete driver;delete example;}

    return 0;
}
//#####################################################################
