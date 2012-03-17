#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <cmath>
#include "ACCURACY_INFO.h"
#include "BOUNDARY_CONDITIONS_BOX.h"
#include "BOUNDARY_CONDITIONS_CIRCLE.h"
#include "BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE.h"
#include "HEADER.h"
#include "TEST_COMMON.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef VECTOR<T,2> TV;
    const int d=TV::m;
    typedef VECTOR<int,d> TV_INT;

    TEST_COMMON<TV> tc;
    tc.Init_1();
    tc.parse_args.Add_String_Argument("-bc_types","sfwf","[swf][swf][swf][swf] source/wall/free for left, bottom, right, top");
    tc.parse_args.Add_Double_Argument("-offset",0,"Offset bottom face");
    tc.parse_args.Add_Double_Argument("-angle",0,"Angle top and bottom faces");
    tc.parse_args.Add_Option_Argument("-gibou_analytic_one","");
    tc.parse_args.Add_Option_Argument("-circle","Circular boundaray");
    tc.parse_args.Add_String_Argument("-error_image","","Image showing error distribution");

    tc.Init_2(argc,argv);
    output_directory=tc.output_directory;

    std::string bc_types=tc.parse_args.Get_String_Value("-bc_types");
    PHYSBAM_ASSERT(bc_types.length()==4 && bc_types.find_first_not_of("swf")==std::string::npos);
    T offset=tc.parse_args.Get_Double_Value("-offset");
    T angle=tc.parse_args.Get_Double_Value("-angle")*(pi/180);
    std::string error_image=tc.parse_args.Get_String_Value("-error_image");
    TV_INT error_image_dim;
    sscanf(error_image.c_str(),"%dx%d",&error_image_dim.x,&error_image_dim.y);

    if(tc.parse_args.Is_Value_Set("-circle"))
        tc.sim.obj.bc=new BOUNDARY_CONDITIONS_CIRCLE<TV>(tc.sim.obj.grid,bc_types[0]=='f'?BOUNDARY_CONDITIONS<TV>::free:bc_types[0]=='w'?BOUNDARY_CONDITIONS<TV>::noslip:BOUNDARY_CONDITIONS<TV>::none);
    else if(tc.parse_args.Is_Value_Set("-gibou_analytic_one")){/*use_accuracy_samples=false;use_extrapolation=false;*/tc.sim.obj.bc=new BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<TV>(tc.sim.obj.grid);}
    else tc.sim.obj.bc=new BOUNDARY_CONDITIONS_BOX<TV>(tc.sim.obj.grid,offset,angle,bc_types);

    tc.Init_3();

    tc.sim.obj.bc->Update_Parameters(tc.sim.param);
    ARRAY<T,FACE_INDEX<d> > u(tc.sim.obj.grid,4),u2(tc.sim.obj.grid,4);
    tc.sim.obj.bc->Initialize_Velocity_Field(u,tc.sim.param.time);

    Dump_Frame<RW>(u,"init");

    tc.sim.obj.ai.Print("INIT",u);

    for(int k=0;k<tc.sim.steps;k++){
        First_Order_Step(tc.sim,u);
//        Second_Order_RE_Step(tc.sim,u,u2);
        Dump_Frame<RW>(u,"step");
        Dump_Error(tc.sim,u,u2);
        tc.sim.obj.bc->Compute_Error(u,tc.sim.param.time);
        if(error_image!="") Dump_Error_Image(tc.sim,u,error_image_dim);}

    tc.sim.obj.ai.Print("END",u);
/*
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(tc.sim.obj.grid,0,GRID<TV>::WHOLE_REGION,-1,0);it.Valid();it.Next()){
        if(it.index.y==1 && it.index.x>1) LOG::cout<<std::endl;
        TV X(tc.sim.obj.grid.Axis_X_Face(it.Full_Index()));
        T e=0;
        if(tc.sim.obj.bc->Inside(X)){
            TV p=tc.sim.obj.bc->Analytic_Velocity(X,tc.sim.param.time);
            e=abs(u(it.Full_Index())-p(it.Axis()));}
        LOG::cout<<it.index.x<<" "<<it.index.y<<" "<<e<<std::endl;}
*/
    return 0;
}
