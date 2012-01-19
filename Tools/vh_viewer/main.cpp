//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "VIEWER_WINDOW.h"
#include <Fl/Fl.H>
#include <Fl/Fl_Window.H>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
using namespace PhysBAM;

int main(int argc,char *argv[])
{
    std::cout<<"Welcome to VH Viewer!"<<std::endl;
    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-offscreen");
    parse_args.Add_String_Argument("-o","","output filename","output filename");
    parse_args.Add_String_Argument("-d","","data directory","data directory");
    parse_args.Add_Vector_3D_Argument("-box_min",VECTOR<double,3>(),"minimum box dimensions","minimum box dimensions");
    parse_args.Add_Vector_3D_Argument("-box_max",VECTOR<double,3>(),"maximum box dimensions","maximum box dimensions");
    parse_args.Add_Integer_Argument("-granularity",100,"granularity","granularity in all directions");
    parse_args.Set_Extra_Arguments(-1,"indices to be removed","indices to be removed");
    parse_args.Parse(argc,argv);

    std::string output_filename="visible_human_levelset.phi";
    std::string data_directory="../../Private_Data/VH_Raw";
    VECTOR<double,3> box_min,box_max;int granularity=100;
    bool offscreen=parse_args.Get_Option_Value("-offscreen");
    if(parse_args.Is_Value_Set("-o")) output_filename=parse_args.Get_String_Value("-o");
    if(parse_args.Is_Value_Set("-d")) data_directory=parse_args.Get_String_Value("-d");
    if(parse_args.Is_Value_Set("-box_min")) box_min=parse_args.Get_Vector_3D_Value("-box_min");
    if(parse_args.Is_Value_Set("-box_max")) box_max=parse_args.Get_Vector_3D_Value("-box_max");
    if(parse_args.Is_Value_Set("-granularity")) granularity=parse_args.Get_Integer_Value("-granularity");

    std::cout<<"Data directory: "<<data_directory<<", output file: "<<output_filename<<std::endl;
    if(offscreen){
        VH_LEVELSET_BUILDER<float> visible_human(data_directory);
        ARRAY<int> tissues;visible_human.All_Tissues_Array(tissues);
        for(int i=0;i<parse_args.Num_Extra_Args();i++) tissues.Remove_Index(atoi(parse_args.Extra_Arg(i).c_str()));
        RANGE<VECTOR<int,3> > box(visible_human.Get_Bounding_Box(tissues));
        if(box_min.x){box.min_corner.x=(int)box_min.x;}if(box_min.y){box.min_corner.y=(int)box_min.y;}if(box_min.z){box.min_corner.z=(int)box_min.z;}
        if(box_max.x){box.max_corner.x=(int)box_max.x;}if(box_max.y){box.max_corner.y=(int)box_max.y;}if(box_max.z){box.max_corner.z=(int)box_max.z;}
        VECTOR<int,3> grid_size=visible_human.Get_Grid_Size(box);
        double scale_factor=(double)granularity/100;VECTOR<int,3> granularity_xyz((int)(grid_size.x*scale_factor),(int)(grid_size.y*scale_factor),(int)(grid_size.z*scale_factor));
        visible_human.Create_Levelset(tissues,box,granularity_xyz,output_filename);
    }
    else{
        VIEWER_WINDOW<float> window;
        Fl::run();
    }
    return 0;
}

