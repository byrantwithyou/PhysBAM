//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Tools/Parsing/STRING_UTILITIES.h>
#include "VH_LEVELSET_BUILDER.h"
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    VH_LEVELSET_BUILDER<float> builder;

    PARAMETER_LIST list;list.Begin_Parse(argc,argv);
    ARRAY<int> tissues;
    STRING_UTILITIES::Parse_Integer_List(list.Get_Parameter("tissues",(std::string)"1"),tissues);
    std::string output_filename=list.Get_Parameter("filename",(std::string)"output.phi");
    builder.data_directory=list.Get_Parameter("data",(std::string)builder.data_directory);
    builder.box_padding_voxels=list.Get_Parameter("box_padding_voxels",(int)0);
    builder.box_padding=list.Get_Parameter("box_padding",(float)0);
    bool dump_images=list.Get_Parameter("dump_images",(bool)false);
    list.End_Parse();

    if(dump_images) builder.Dump_Images();
    else builder.Create_Levelset(tissues,output_filename);

    return 0;
}

