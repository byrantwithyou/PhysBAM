#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include "ACCURACY_INFO.h"
#include "BOUNDARY_CONDITIONS_SEGMENT.h"
#include "HEADER.h"
#include "TEST_COMMON.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef VECTOR<T,1> TV;
    const int d=TV::m;
    typedef VECTOR<int,d> TV_INT;

    T x0=0,x1=1;
    TEST_COMMON<TV> tc;
    tc.Init_1();
    tc.parse_args.Add("-x0",&x0,"value","left endpoint (<= 0)");
    tc.parse_args.Add("-x1",&x1,"value","right endpoint (>= 1)");

    tc.Init_2(argc,argv);
    output_directory=tc.output_directory;

    tc.sim.obj.bc=new BOUNDARY_CONDITIONS_SEGMENT<TV>(tc.sim.obj.grid,x0,x1);

    tc.Init_3();

    tc.sim.obj.bc->Update_Parameters(tc.sim.param);
    ARRAY<T,FACE_INDEX<d> > u(tc.sim.obj.grid,4),u2(tc.sim.obj.grid,4);
    tc.sim.obj.bc->Initialize_Velocity_Field(u,tc.sim.param.time);

    Dump_Frame<RW>(u,"init");

    tc.sim.obj.ai.Print("INIT",u);

    for(int k=0;k<tc.sim.steps;k++){
        Second_Order_RE_Step(tc.sim,u,u2);
        Dump_Error(tc.sim,u,u2);}

    tc.sim.obj.ai.Print("END",u);
    tc.sim.obj.bc->Compute_Error(u,tc.sim.param.time);

    return 0;
}
