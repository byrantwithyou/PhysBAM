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
    parse_args.Add_Integer_Argument("-scale",50,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-subscale",2,"fine/coarse scale grid resolution ratio");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
    parse_args.Add_Option_Argument("-binary","use binary refinement");
    parse_args.Add_Double_Argument("-alpha",1,"interpolation parameter");    

    parse_args.Parse();
    parse_args.Print_Arguments();

    SMOKE_TESTS<TV>* example=new SMOKE_TESTS<TV>(stream_type,parse_args);
    INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
