#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
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

    TEST_COMMON<TV> tc;
    tc.Init_1();
    tc.parse_args.Add_Double_Argument("-x0",(T)0,"left endpoint (<= 0)");
    tc.parse_args.Add_Double_Argument("-x1",(T)1,"right endpoint (>= 1)");

    tc.Init_2(argc,argv);
    output_directory=tc.output_directory;

    T x0=tc.parse_args.Get_Double_Value("-x0");
    T x1=tc.parse_args.Get_Double_Value("-x1");

    tc.sim.obj.bc=new BOUNDARY_CONDITIONS_SEGMENT<TV>(tc.sim.obj.grid,x0,x1);

    tc.Init_3();

    tc.sim.obj.bc->Update_Parameters(tc.sim.param);
    ARRAY<T,FACE_INDEX<d> > u(tc.sim.obj.grid,4),u2(tc.sim.obj.grid,4);
    tc.sim.obj.bc->Initialize_Velocity_Field(u,tc.sim.param.time);

    Dump_Frame<RW>(u,"init");

    tc.sim.obj.ai.Print("INIT",u);

    for(int k=1;k<=tc.sim.steps;k++){
        Second_Order_RE_Step(tc.sim,u,u2);
        Dump_Error(tc.sim,u,u2);}

    tc.sim.obj.ai.Print("END",u);
    tc.sim.obj.bc->Compute_Error(u,tc.sim.param.time);

    return 0;
}
