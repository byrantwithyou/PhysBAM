//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Images/DEPTH_BUFFERING.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Dynamics/Read_Write/EPS_FILE_GEOMETRY.h>

using namespace PhysBAM;

typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> TV_INT;
typedef VECTOR<T,2> TV2;
typedef VECTOR<int,2> TV_INT2;
typedef DEPTH_BUFFERING<T> DB;

int main()
{
    DEPTH_BUFFERING<T> db;
    db.Add_Element(TV(0,0,0),TV(1,0,0),TV(0,1,1),0);
    db.Add_Element(TV(0,0,1),TV(1,0,1),TV(0,1,0),1);
    ARRAY<DISPLAY_PRIMITIVE_ORDERING<T> >& dps=db.Process_Primitives();
    EPS_FILE_GEOMETRY<T> f("out.eps");

    for(int i=0;i<dps.m;i++){
        TRIANGLE_2D<T> triangle(DB::Project(dps(i).vertices(0)),DB::Project(dps(i).vertices(1)),DB::Project(dps(i).vertices(2)));
        switch(dps(i).style){
            case 0:
                f.Fill_Object(triangle,"0 0 1");
                break;
            case 1:
                // f.Draw_Object(triangle);
                break;
            default: break;}}
    return 0;
}
