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

typedef double RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,2> TV_INT;

const int w=3;

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    RANDOM_NUMBERS<T> random;
    random.Set_Seed(0);

    MPLE_DRIVER<TV,w> mple;
    mple.grid.Initialize(TV_INT(50+1,50+1),RANGE<TV>(TV(-1,-1),TV(1,1)));
    std::string output_directory="output";
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),mple.grid,output_directory);

    const int number_of_points=2; 
    const RANGE<TV> block(TV(-.1,-.1),TV(.1,.1));
    
    for(int i=0;i<number_of_points;i++){
        MPLE_POINT<TV,w> point;
        point.X=random.Get_Uniform_Vector(block);
        mple.points.Append(point);}

    mple.Run();

    return 0;
}
