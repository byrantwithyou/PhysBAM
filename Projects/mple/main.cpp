//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
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

    RANDOM_NUMBERS<T> random;
    random.Set_Seed(0);

    MPLE_DRIVER<TV,w> mple;
    const RANGE<TV> domain(TV(-.5,-.5),TV(.5,.5));
    const TV_INT counts(TV_INT(50,50));
    mple.grid.Initialize(counts+1,domain);

    std::string output_directory="output";
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),mple.grid,output_directory);

    const int number_of_points=1000; 
    const RANGE<TV> block(TV(-.3,-.3),TV(.3,.3));
    
    for(int i=0;i<number_of_points;i++){
        MPLE_POINT<TV,w> point;
        point.X=random.Get_Uniform_Vector(block);
        mple.points.Append(point);}

    mple.Run();

    LOG::Finish_Logging();
    return 0;
}
