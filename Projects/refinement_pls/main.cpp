#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Dynamics/PLS_EXAMPLE.h>
#include <PhysBAM_Dynamics/PLS_REFINEMENT_DRIVER.h>
#include "Water_Tests/WATER_TESTS.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    LOG::Initialize_Logging(false,false,1<<30,true);
    
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
//    typedef VECTOR<T,2> TV;
    typedef VECTOR<T,3> TV;

    MPI_WORLD mpi_world(argc,argv);

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Integer_Argument("-restart",0,"restart frame");
    parse_args.Add_String_Argument("-split","","split restart data");
    parse_args.Add_Integer_Argument("-scale",128,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-subscale",2,"fine/coarse scale grid resolution ratio");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
    parse_args.Add_Integer_Argument("-threads",0,"number of threads");
    parse_args.Add_Integer_Argument("-buffer",1,"surface buffer");
    parse_args.Add_Double_Argument("-cfl",.9,"CFL");
    parse_args.Add_Double_Argument("-x",1,"ratio multiplier");
    parse_args.Add_Double_Argument("-y",1,"ratio multiplier");
    parse_args.Add_Double_Argument("-z",1,"ratio multiplier");
    parse_args.Add_String_Argument("-output","output","output directory");
    parse_args.Add_Integer_Argument("-binary",0,"total number of binary adaptive refinement levels");
    parse_args.Add_Double_Argument("-alpha",1,"interpolation parameter");
    parse_args.Add_Integer_Argument("-last_frame",100,"number of frames simulated");
    parse_args.Add_Option_Argument("-3d","3D examples");
    parse_args.Add_Option_Argument("-cubic","cubic interpolation");
    parse_args.Add_Option_Argument("-nosurface","don't solve the surface as one solve");
    parse_args.Add_Option_Argument("-write_debug","write debug data");
    parse_args.Add_Integer_Argument("-scheme",1,"scheme type for binary");

    parse_args.Parse();
    parse_args.Print_Arguments();
    
    WATER_TESTS<TV>* example=new WATER_TESTS<TV>(stream_type,parse_args);

    if(mpi_world.initialized){
        example->coarse_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->coarse_mac_grid,3);
        example->fine_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->fine_mac_grid,3);
        if(example->fine_mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->fine_mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}
    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    PLS_REFINEMENT_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
