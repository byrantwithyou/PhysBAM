//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include "MPLE_DRIVER.h"

using namespace PhysBAM;

const int w=4;
const int d=2;

typedef float RW;
typedef double T;
typedef VECTOR<T,d> TV;
typedef VECTOR<int,d> TV_INT;

int main(int argc,char* argv[])
{
    LOG::Initialize_Logging();
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse(true);

    bool use_test=false;
    int test_number=-1;
    bool use_output_directory=false;
    std::string output_directory="output";
    MPLE_DRIVER<TV,w> mple;
    RANDOM_NUMBERS<T> random;
    random.Set_Seed(0);
    int resolution=32;
    int number_of_points=250;
    
    parse_args.Add("-test",&test_number,&use_test,"test","test number");
    parse_args.Add("-o",&output_directory,&use_output_directory,"o","output directory");
    parse_args.Add("-spread",&mple.spread,"spread","Interface spread int cells");
    parse_args.Add("-contour_value",&mple.contour_value,"contour value","Contour value");
    parse_args.Add("-mu",&mple.mu,"mu","Source power");
    parse_args.Add("-frames",&mple.frames,"frames","Number of frames");
    parse_args.Add("-cfl",&mple.cfl,"cfl","CFL number");
    parse_args.Add("-rescale",&mple.rescale,"rescale","Rescale and clamp source");
    parse_args.Add("-identity",&mple.identity,"identity","Identity scale");
    parse_args.Add("-resolution",&resolution,"resolution","Grid resolution");
    parse_args.Add("-points",&number_of_points,"points","Number of sample points");
    parse_args.Parse(true);
    parse_args.Parse();

    if(!use_test){LOG::cerr<<"Test number is required."<<std::endl;exit(-1);}
    if(!use_output_directory) output_directory="output";
    mple.output_directory=output_directory;

    switch(test_number){
        case 1:{
            const RANGE<TV> domain(TV(-.75,-.75),TV(.75,.75));
            const TV_INT counts(TV_INT(resolution,resolution));
            mple.grid.Initialize(counts+1,domain);
            const RANGE<TV> block(TV(-.3,-.3),TV(.3,.3));
            const ROTATION<TV> rotation=ROTATION<TV>::From_Angle(45);
            for(int i=0;i<number_of_points;i++){
                MPLE_POINT<TV,w> point;
                point.X=rotation.Rotate(random.Get_Uniform_Vector(block));
                mple.points.Append(point);}}
            break;
        default:{LOG::cerr<<"Unknown test number"<<std::endl;exit(-1);}}

    LOG::cout<<"### PARAMETERS ###"<<std::endl;
    LOG::cout<<"cfl "<<mple.cfl<<std::endl;
    LOG::cout<<"frames "<<mple.frames<<std::endl;
    LOG::cout<<"mu "<<mple.mu<<std::endl;
        
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),mple.grid,output_directory);
    mple.Run();

    LOG::Finish_Logging();
    return 0;
}
