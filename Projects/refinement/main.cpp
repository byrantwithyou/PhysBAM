#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_REFINEMENT_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_REFINEMENT_EXAMPLE.h>
#include "Smoke_Tests/SMOKE_TESTS.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    LOG::Initialize_Logging(false,false,1<<31,true);

    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    MPI_WORLD mpi_world(argc,argv);

    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-restart",0,"restart frame");
    parse_args.Add_String_Argument("-split","","split restart data");
    parse_args.Add_Integer_Argument("-scale",64,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-subscale",4,"fine/coarse scale grid resolution ratio");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
    parse_args.Add_Integer_Argument("-threads",0,"number of threads");
    parse_args.Add_Integer_Argument("-x",1,"ratio multiplier");
    parse_args.Add_Integer_Argument("-y",2,"ratio multiplier");
    parse_args.Add_Integer_Argument("-z",1,"ratio multiplier");
    parse_args.Add_Option_Argument("-coarse_forces","Use coarse forces");
    parse_args.Add_Double_Argument("-cfl",.9,"CFL");
    parse_args.Add_String_Argument("-output","output","output directory");
    parse_args.Add_Integer_Argument("-binary",0,"total number of binary adaptive refinement levels");
    parse_args.Add_Double_Argument("-alpha",1,"interpolation parameter");
    parse_args.Add_Double_Argument("-vc",.1,"vorticity confinement");
    parse_args.Add_Double_Argument("-buoyancy",1,"buoyancy constant");
    parse_args.Add_Integer_Argument("-last_frame",200,"number of frames simulated");
    parse_args.Add_Integer_Argument("-scheme",2,"scheme type for binary");
    parse_args.Add_Integer_Argument("-test_number",1,"test number");
    parse_args.Add_Double_Argument("-kol",0,"kolmogorov");
    parse_args.Add_Integer_Argument("-short_kol",0,"turn off kolmogorov after the specified frame. 0 means don't turn off.");
    parse_args.Add_Double_Argument("-source_radius",0.05,"radius of source");
    parse_args.Add_Option_Argument("-3d","do 3d solver");
    
    parse_args.Parse(argc,argv);
    parse_args.Print_Arguments(argc,argv);

    if(!parse_args.Is_Value_Set("-3d")){
        typedef VECTOR<T,2> TV;
        typedef VECTOR<int,TV::dimension> TV_INT;

        SMOKE_TESTS<TV>* example=new SMOKE_TESTS<TV>(stream_type,parse_args);

        if(mpi_world.initialized){
            example->coarse_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->coarse_mac_grid,3);
            example->fine_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->fine_mac_grid,3);
            if(example->fine_mpi_grid->Number_Of_Processors()>1){
                example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->fine_mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}}
        FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
        
        INCOMPRESSIBLE_REFINEMENT_DRIVER<TV> driver(*example);
        driver.Execute_Main_Program();}
    else{
        typedef VECTOR<T,3> TV;
        typedef VECTOR<int,TV::dimension> TV_INT;

        SMOKE_TESTS<TV>* example=new SMOKE_TESTS<TV>(stream_type,parse_args);

        if(mpi_world.initialized){
            example->coarse_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->coarse_mac_grid,3);
            example->fine_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->fine_mac_grid,3);
            if(example->fine_mpi_grid->Number_Of_Processors()>1){
                example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->fine_mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}}
        FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
        
        INCOMPRESSIBLE_REFINEMENT_DRIVER<TV> driver(*example);
        driver.Execute_Main_Program();}

    return 0;
}
