#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_ADAPTIVE_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_ADAPTIVE_EXAMPLE.h>
#include "Smoke_Tests/SMOKE_TESTS.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;
//    typedef VECTOR<T,3> TV;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Print_Arguments();

    SMOKE_TESTS<TV>* example=new SMOKE_TESTS<TV>(stream_type,parse_args);
    INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
